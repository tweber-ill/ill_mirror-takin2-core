/**
 * Takin (TAS tool)
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2014
 * @license GPLv2
 */

#include "taz.h"
#include "tlibs/string/spec_char.h"
#include "tlibs/log/log.h"
#include "tlibs/math/rand.h"
#include "tlibs/log/debug.h"
#include "tlibs/time/chrono.h"
#include "libs/globals.h"
#include "dialogs/NetCacheDlg.h"

#include "tools/monteconvo/ConvoDlg.h"
#include "tools/monteconvo/monteconvo_cli.h"
#include "tools/convofit/convofit_cli.h"
#include "tools/scanviewer/scanviewer.h"

#include "tlibs/version.h"
#include "libs/version.h"
#include "libcrystal/version.h"

#include <system_error>
#include <boost/version.hpp>
#include <boost/system/system_error.hpp>
#include <boost/asio/io_service.hpp>
#include <boost/asio/signal_set.hpp>
#include <boost/scope_exit.hpp>
#include <boost/program_options.hpp>

#include <locale>
#include <clocale>
#include <memory>

#include <QMetaType>
#include <QTextCursor>
#include <QDir>
#include <QMessageBox>
#include <QSplashScreen>
#include <QStyleFactory>

#include <unistd.h>


namespace chr = std::chrono;
namespace asio = boost::asio;
namespace sys = boost::system;
namespace opts = boost::program_options;


// ----------------------------------------------------------------------------
// hacks
#ifdef NON_STANDALONE_MINUIT
	//#include <root/TError.h>
	void SetErrorHandler(void (*)(int, bool, const char*, const char*));
#else
	// dummy handler
	void SetErrorHandler(void (*)(int, bool, const char*, const char*)) {}
#endif
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// logging
static bool add_logfile(std::ofstream* postrLog, bool bAdd=1)
{
	if(!postrLog || !postrLog->is_open())
	{
		tl::log_err("Cannot open log file.");
		return 0;
	}

	for(tl::Log* plog : { &tl::log_info, &tl::log_warn, &tl::log_err, &tl::log_crit, &tl::log_debug })
	{
		if(bAdd)
			plog->AddOstr(postrLog, 0, 0);
		else
			plog->RemoveOstr(postrLog);
	}

	if(!bAdd) postrLog->operator<<(std::endl);
	return 1;
}

template<class SysErr = std::system_error>
static inline void sys_err(const SysErr& err)
{
	tl::log_crit("System error: ", err.what(),
		", category: ", err.code().category().name(),
		", value: ", err.code().value(), ".");
	tl::log_backtrace();
}


static void show_splash_msg(QApplication *pApp, QSplashScreen *pSplash, const std::string &strMsg)
{
	if(!pApp || !pSplash)
		return;

	QColor colSplash(0xff, 0xcc, 0x00);
	pSplash->showMessage(strMsg.c_str(), Qt::AlignCenter, colSplash);
	pApp->processEvents();
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
class TakAppl : public QApplication
{
protected:
	std::shared_ptr<TazDlg> m_pTakDlg;
	std::string m_strToLoad;

protected:

public:
	// need a reference to argc, because the QApplication constructor
	// would otherwise take the temporary stack variable
	// see: https://doc.qt.io/qt-5/qapplication.html#QApplication
	TakAppl(/*const*/ int& argc, /*const*/ char** argv) : QApplication(argc, argv) {}
	virtual ~TakAppl() {}

	void SetTakDlg(std::shared_ptr<TazDlg> pDlg) { m_pTakDlg = pDlg; }


	void DoPendingRequests()
	{
		if(!m_pTakDlg) return;

		// user clicked on an associated file to load?
		if(m_strToLoad != "")
		{
			m_pTakDlg->Load(m_strToLoad.c_str());
		}
	}


	virtual bool event(QEvent *pEvt) override
	{
		if(pEvt->type() == QEvent::FileOpen)
		{
			m_strToLoad = ((QFileOpenEvent*)pEvt)->file().toStdString();
			return true;
		}

		return QApplication::event(pEvt);
	}


	virtual bool notify(QObject* obj, QEvent* evt) override
	{
		try
		{
			return QApplication::notify(obj, evt);
		}
		catch(const std::exception& ex)
		{
			tl::log_err(ex.what());
		}

		return false;
	}
};

// ----------------------------------------------------------------------------



/**
 * entry point for the main GUI
 */
int main(int argc, char** argv)
{
	try
	{
		pid_t pidMain = getpid();
		std::ios_base::sync_with_stdio(0);


#ifdef NO_TERM_CMDS
		tl::Log::SetUseTermCmds(0);
#endif

		// install exit signal handlers
		asio::io_service ioSrv;
		asio::signal_set sigInt(ioSrv, SIGABRT, SIGTERM, SIGINT);
		sigInt.async_wait([&ioSrv, pidMain](const sys::error_code& err, int iSig)
		{
			pid_t pid = getpid();
			if(pid != pidMain)
			{
				tl::log_warn("Child process exit requested via signal.", pid, ".");
				ioSrv.stop();
				exit(-1);
			}

			tl::log_warn("Hard exit requested via signal ", iSig, ". This may cause a fault.");
			if(err) tl::log_err("Error: ", err.message(), ", error category: ", err.category().name(), ".");
			ioSrv.stop();
#ifdef SIGKILL
			std::raise(SIGKILL);
#endif
			exit(-1);
		});
		std::thread thSig([&ioSrv]() { ioSrv.run(); });
		BOOST_SCOPE_EXIT(&ioSrv, &thSig)
		{
			ioSrv.stop();
			thSig.join();
		}
		BOOST_SCOPE_EXIT_END
		// --------------------------------------------------------------------


		// only for non-standalone minuit
		SetErrorHandler([](int, bool, const char*, const char* pcMsg) { tl::log_err(pcMsg); });



		// --------------------------------------------------------------------
		// get program options
		std::vector<std::string> vecTazFiles;
		bool bShowHelp = 0;
		bool bStartScanviewer = 0;
		bool bStartMonteconvo = 0;
		bool bStartMonteconvoCLI = 0;
		bool bStartConvofit = 0;
		bool bStartTakinMain = 1, bStartGui = 1;

		opts::options_description args("Takin options");
		args.add(boost::shared_ptr<opts::option_description>(
			new opts::option_description("taz-file",
			opts::value<decltype(vecTazFiles)>(&vecTazFiles),
			"loads a Takin session file")));
		args.add(boost::shared_ptr<opts::option_description>(
			new opts::option_description("help",
			opts::bool_switch(&bShowHelp),
			"shows the program options")));
		args.add(boost::shared_ptr<opts::option_description>(
			new opts::option_description("scans",
			opts::bool_switch(&bStartScanviewer),
			"directly runs the scan viewer tool")));
		args.add(boost::shared_ptr<opts::option_description>(
			new opts::option_description("convo",
			opts::bool_switch(&bStartMonteconvo),
			"directly runs the convolution simulator/fitter")));
		args.add(boost::shared_ptr<opts::option_description>(
			new opts::option_description("convosim",
			opts::bool_switch(&bStartMonteconvoCLI),
			"runs the convolution simulation command-line tool")));
		args.add(boost::shared_ptr<opts::option_description>(
			new opts::option_description("convofit",
			opts::bool_switch(&bStartConvofit),
			"runs the convolution fitting command-line tool")));

		// positional args
		opts::positional_options_description args_pos;
		args_pos.add("taz-file", -1);

		opts::basic_command_line_parser<char> clparser(argc, argv);
		clparser.options(args);
		clparser.positional(args_pos);
		opts::basic_parsed_options<char> parsedopts = clparser.run();

		opts::variables_map opts_map;
		opts::store(parsedopts, opts_map);
		opts::notify(opts_map);

		if(bStartMonteconvo || bStartScanviewer || bStartMonteconvoCLI || bStartConvofit)
			bStartTakinMain = 0;
		if(bStartMonteconvoCLI || bStartConvofit)
			bStartGui = 0;

		if(bShowHelp)
		{
			std::cout << args << std::endl;
			return 0;
		}
		// --------------------------------------------------------------------



		std::string strLog = QDir::tempPath().toStdString();
		strLog += "/takin.log";
		std::ofstream ofstrLog(strLog, std::ios_base::out|std::ios_base::app);
		if(add_logfile(&ofstrLog, 1))
			tl::log_debug("Logging to file \"", strLog, "\".");


		tl::log_info("--------------------------------------------------------------------------------");
		tl::log_info("This is Takin version " TAKIN_VER ".");
		tl::log_info("Author: Tobias Weber <tweber@ill.fr>, 2014 - 2020.");
		tl::log_info("Licensed under GPLv2, see the about dialog.");
		tl::log_info("--------------------------------------------------------------------------------");

		tl::log_debug("Using ", sizeof(t_real_glob)*8, " bit ", tl::get_typename<t_real_glob>(), "s as internal data type.");


		std::unique_ptr<QCoreApplication> app;
		TakAppl *app_gui = nullptr;

		if(bStartGui)
			app.reset(app_gui = new TakAppl(argc, argv));
		else
			app.reset(new QCoreApplication(argc, argv));

		app->setApplicationName("Takin");
		app->setApplicationVersion(TAKIN_VER);


		// qt needs to be able to copy these structs when emitting signals from a different thread
		qRegisterMetaType<TriangleOptions>("TriangleOptions");
		qRegisterMetaType<CrystalOptions>("CrystalOptions");
		qRegisterMetaType<std::string>("std::string");
		qRegisterMetaType<CacheVal>("CacheVal");
		qRegisterMetaType<QTextCursor>("QTextCursor");


		// locale
		std::setlocale(LC_ALL, "C");
		std::locale::global(std::locale::classic());
		QLocale::setDefault(QLocale::English);


		tl::init_rand();

		g_strHome = QDir::homePath().toStdString() + "/.takin";
		g_strApp = QCoreApplication::applicationDirPath().toStdString();
		tl::log_info("Program path: ", g_strApp);
		tl::log_info("Home path: ", g_strHome);

		add_resource_path(g_strHome, 0);
		add_resource_path(g_strApp);
		add_resource_path(g_strApp + "/..");
		add_resource_path(g_strApp + "/resources");
		add_resource_path(g_strApp + "/Resources");
		add_resource_path(g_strApp + "/../resources");
		add_resource_path(g_strApp + "/../Resources");

		QCoreApplication::addLibraryPath((g_strApp + "/../lib/plugins").c_str());
		QCoreApplication::addLibraryPath((g_strApp + "/lib/plugins").c_str());


		// run command-line tools
		if(bStartMonteconvoCLI)
			return monteconvo_main(argc, argv);
		else if(bStartConvofit)
			return convofit_main(argc, argv);


		// ------------------------------------------------------------
		// GUI stuff

		#if defined Q_WS_X11 && !defined NO_3D
			//XInitThreads();
			QCoreApplication::setAttribute(Qt::AA_X11InitThreads, true);
			QGL::setPreferredPaintEngine(QPaintEngine::OpenGL);
		#endif

		QSettings settings("tobis_stuff", "takin");


		// set user-selected GUI style
		if(settings.contains("main/gui_style_value"))
		{
			QString strStyle = settings.value("main/gui_style_value", "").toString();
			QStyle *pStyle = QStyleFactory::create(strStyle);
			if(pStyle)
				QApplication::setStyle(pStyle);
			else
				tl::log_err("Style \"", strStyle.toStdString(), "\" was not found.");
		}


		// ------------------------------------------------------------
		// splash screen
		std::unique_ptr<QSplashScreen> pSplash;

		if(bStartTakinMain)
		{
			QPixmap pixSplash = load_pixmap("res/icons/takin.svg");
			if(!pixSplash.isNull())
				pixSplash = pixSplash.scaled(pixSplash.size().width()*0.55, pixSplash.size().height()*0.55);
			pSplash.reset(new QSplashScreen{pixSplash});
		}

		if(pSplash)
		{
			QFont fontSplash = pSplash->font();
			fontSplash.setPixelSize(14);
			fontSplash.setBold(1);
			pSplash->setFont(fontSplash);
			pSplash->show();
		}

		const std::string strStarting = "Starting up Takin version " TAKIN_VER ".";
		show_splash_msg(app_gui, pSplash.get(), strStarting);


		// ------------------------------------------------------------
		// tlibs & libcrystal version checks
		tl::log_info("Using tlibs version ", tl::get_tlibs_version(), ".");
		if(!tl::check_tlibs_version(TLIBS_VERSION))
		{
			tl::log_crit("Version mismatch in tLibs. Please recompile.");
			tl::log_crit("tLibs versions: library: ", tl::get_tlibs_version(),
				", headers: ", TLIBS_VERSION, ".");

			QMessageBox::critical(0, "Takin - Error", "Broken build: Mismatch in tlibs version.");
			return -1;
		}

		tl::log_info("Using libcrystal version ", xtl::get_libcrystal_version(), ".");
		if(!xtl::check_libcrystal_version(LIBCRYSTAL_VERSION))
		{
			tl::log_crit("Version mismatch in LibCrystal. Please recompile.");
			tl::log_crit("LibCrystal versions: library: ", xtl::get_libcrystal_version(),
				", headers: ", LIBCRYSTAL_VERSION, ".");

			QMessageBox::critical(0, "Takin - Error", "Broken build: Mismatch in libcrystal version.");
			return -1;
		}


		#define TAKIN_CHECK " Please check if Takin is correctly installed and the current working directory is set to the Takin main directory."
		show_splash_msg(app_gui, pSplash.get(), strStarting + "\nChecking resources ...");

		// check tables
		g_bHasElements = (find_resource("res/data/elements.xml") != "");
		g_bHasScatlens = (find_resource("res/data/scatlens.xml") != "");
		g_bHasFormfacts = (find_resource("res/data/ffacts.xml") != "");
		g_bHasMagFormfacts = (find_resource("res/data/magffacts.xml") != "");
		g_bHasSpaceGroups = (find_resource("res/data/sgroups.xml") != "");

		if(!g_bHasScatlens)
		{
			const char* pcErr = "Scattering length table could not be found." TAKIN_CHECK;
			tl::log_err(pcErr);

			QMessageBox::critical(0, "Takin - Error", pcErr);
			return -1;
		}
		if(!g_bHasFormfacts)
		{
			const char* pcErr = "Atomic form factor coefficient table could not be found." TAKIN_CHECK;
			tl::log_err(pcErr);

			QMessageBox::critical(0, "Takin - Error", pcErr);
			return -1;
		}
		if(!g_bHasMagFormfacts)
		{
			const char* pcErr = "Magnetic form factor coefficient table could not be found." TAKIN_CHECK;
			tl::log_warn(pcErr);

			//QMessageBox::warning(0, "Takin - Warning", pcErr);
			//return -1;
		}
		if(!g_bHasSpaceGroups)
		{
			const char* pcErr = "Space group table could not be found!" TAKIN_CHECK;
			tl::log_err(pcErr);

			QMessageBox::critical(0, "Takin - Error", pcErr);
			return -1;
		}

		tl::init_spec_chars();


		// check if icons are available
		if(find_resource("res/icons/document-new.svg") == "" ||
			find_resource("res/icons/takin.svg") == "")
		{
			const char* pcErr = "Takin resources could not be found!" TAKIN_CHECK;
			tl::log_err(pcErr);

			QMessageBox::critical(0, "Takin - Error", pcErr);
			return -1;
		}


		// ------------------------------------------------------------
		std::shared_ptr<TazDlg> pTakDlg;
		std::shared_ptr<ConvoDlg> pConvoDlg;
		std::shared_ptr<ScanViewerDlg> pScanViewerDlg;

		if(bStartTakinMain)
		{
			{
	#ifdef IS_EXPERIMENTAL_BUILD
				int iPrevDaysSinceEpoch = 0;
				if(settings.contains("debug/last_warned"))
					iPrevDaysSinceEpoch = settings.value("debug/last_warned").toInt();
				int iDaysSinceEpoch = tl::epoch_dur<tl::t_dur_days<int>>().count();

				std::string strExp = "This " BOOST_PLATFORM " version of Takin is still experimental, "
					"does not include all features and may show unexpected behaviour. Please report "
					"bugs to tobias.weber@tum.de. Thanks!";
				tl::log_warn(strExp);
				//tl::log_debug("Days since last warning: ", iDaysSinceEpoch-iPrevDaysSinceEpoch, ".");

				// show warning message box every 5 days
				if(iDaysSinceEpoch - iPrevDaysSinceEpoch >= 5)
				{
					QMessageBox::warning(0, "Takin", strExp.c_str());
					settings.setValue("debug/last_warned", iDaysSinceEpoch);
				}
	#endif

			/*// Warnings due to version changes
				if(settings.value("debug/last_warning_shown", 0).toInt() < 1)
				{
					QMessageBox::warning(0, "Takin", "Please beware that in this version "
						"the correction factors for the resolution convolution have been re-worked. "
						"Any global scaling factors will have changed.");
					settings.setValue("debug/last_warning_shown", 1);
				}*/
			}


			show_splash_msg(app_gui, pSplash.get(), strStarting + "\nLoading 1/2 ...");
			pTakDlg.reset(new TazDlg{nullptr, strLog});
			app_gui->SetTakDlg(pTakDlg);
			show_splash_msg(app_gui, pSplash.get(), strStarting + "\nLoading 2/2 ...");
			app_gui->DoPendingRequests();

			if(pSplash) pSplash->finish(pTakDlg.get());
			if(vecTazFiles.size() >= 1)
			{
				tl::log_info("Loading \"", vecTazFiles[0], "\"...");
				pTakDlg->Load(vecTazFiles[0].c_str());
			}
			pTakDlg->show();
		}
		if(bStartMonteconvo)
		{
			pConvoDlg.reset(new ConvoDlg{nullptr, &settings});
			pConvoDlg->setWindowFlags(Qt::Window);

			// load a given convolution session file
			if(vecTazFiles.size() >= 1)
			{
				tl::log_info("Loading \"", vecTazFiles[0], "\"...");

				const std::string strXmlRoot("taz/");
				tl::Prop<std::string> xml;
				if(xml.Load(vecTazFiles[0], tl::PropType::XML))
					pConvoDlg->Load(xml, strXmlRoot);
				else
					QMessageBox::critical(nullptr, "Error", "Could not load convolution file.");
			}

			pConvoDlg->show();
		}
		if(bStartScanviewer)
		{
			pScanViewerDlg.reset(new ScanViewerDlg{nullptr});
			pScanViewerDlg->setWindowFlags(Qt::Window);

			pScanViewerDlg->show();
		}

		int iRet = app_gui->exec();

		// ------------------------------------------------------------

		tl::deinit_spec_chars();
		tl::log_info("Shutting down Takin.");
		add_logfile(&ofstrLog, 0);

		return iRet;
	}
	catch(const std::system_error& err) { sys_err(err); }
	catch(const boost::system::system_error& err) { sys_err(err); }
	catch(const std::exception& ex)
	{
		tl::log_crit("Exception: ", ex.what());
		tl::log_backtrace();
	}

	return -1;
}
