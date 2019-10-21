#
# implementation of the eckold-sobolev algo
#
# @author Tobias Weber <tweber@ill.fr>
# @date feb-2015, oct-2019
# @license GPLv2
#
# @desc for algorithm: [eck14] G. Eckold and O. Sobolev, NIM A 752, pp. 54-64 (2014)
# @desc for alternate R0 normalisation: [mit84] P. W. Mitchell, R. A. Cowley and S. A. Higgins, Acta Cryst. Sec A, 40(2), 152-160 (1984)
#

import numpy as np
import numpy.linalg as la


#
# constants
# see: https://code.ill.fr/scientific-software/takin/mag-core/blob/master/tools/tascalc/tascalc.py
#
ksq2E = 2.072124836832
sig2fwhm = 2.*np.sqrt(2.*np.log(2.))


#
# z rotation matrix
#
def rotation_matrix_3d_z(angle):
    if angle == 0.:
        s = 0.
        c = 1.
    else:
        s = np.sin(angle)
        c = np.cos(angle)

    return np.array([
        [c, -s, 0],
        [s,  c, 0],
        [0,  0, 1]
    ])



#
# projects along one axis of the quadric
# see: https://code.ill.fr/scientific-software/takin/mag-core/blob/master/tools/tascalc/cov.py
#
def quadric_proj(_E, idx):
    E = np.delete(np.delete(_E, idx, axis=0), idx, axis=1)
    if np.abs(_E[idx, idx]) < 1e-8:
        return E

    v = 0.5 * (_E[idx,:] + _E[:,idx])
    vv = np.outer(v, v) / _E[idx, idx]
    vv = np.delete(np.delete(vv, idx, axis=0), idx, axis=1)

    return E - vv


#
# Project linear part of the quadric
# (see [eck14], equ. 57)
#
def quadric_proj_vec(vec, _E, idx):
    _col = _E[:,idx]
    col = np.delete(_col, idx, axis=0)
    if np.abs(_col[idx]) < 1e-8:
        return col

    v = np.delete(vec, idx, axis=0)
    v = v - col*vec[idx]/_col[idx]

    return v



def k2lam(k):
    return 2.*np.pi / k



#
# mono (and ana) resolution calculation
#
def get_mono_vals(src_w, src_h, mono_w, mono_h,
	dist_src_mono, dist_mono_sample,
	ki, thetam,
	coll_h_pre_mono, coll_h_pre_sample,
	coll_v_pre_mono, coll_v_pre_sample,
	mono_mosaic, mono_mosaic_v,
	inv_mono_curvh, inv_mono_curvv,
	pos_x , pos_y, pos_z,
	refl):

    # A matrix: formula 26 in [eck14]
    A = np.identity(3)

    A_t0 = 1. / mono_mosaic
    A_tx = inv_mono_curvh*dist_mono_sample / np.abs(np.sin(thetam))
    A_t1 = A_t0*A_tx

    A[0,0] = 0.5*sig2fwhm**2. / ki**2. * np.tan(thetam)**2. * \
        ( (2./coll_h_pre_mono)**2. + (2*dist_src_mono/src_w)**2. + A_t0*A_t0 )

    A[0,1] = A[1,0] = 0.5*sig2fwhm**2. / ki**2. * np.tan(thetam) * \
        ( + 2.*(1./coll_h_pre_mono)**2. + 2.*dist_src_mono*(dist_src_mono-dist_mono_sample)/src_w**2. + \
            A_t0**2. - A_t0*A_t1)

    A[1,1] = 0.5*sig2fwhm**2. / ki**2. * \
    ( (1./coll_h_pre_mono)**2. + (1./coll_h_pre_sample)**2. \
        + ((dist_src_mono-dist_mono_sample)/src_w)**2. \
        + (dist_mono_sample/(mono_w*np.abs(np.sin(thetam))))**2. \
        + A_t0*A_t0 - 2.*A_t0*A_t1 + A_t1*A_t1)



	# Av matrix: formula 38 in [eck14]
	# some typos in paper leading to the (false) result of a better Qz resolution when focusing
	# => trying to match terms in Av with corresponding terms in A
	# corresponding pre-mono terms commented out in Av, as they are not considered there
    Av = np.identity(2)

    Av_t0 = 0.5 / (mono_mosaic_v*np.abs(np.sin(thetam)))
    Av_t1 = inv_mono_curvv*dist_mono_sample / mono_mosaic_v

    Av[0,0] = 0.5*sig2fwhm**2. / ki**2. * \
        ( (1./coll_v_pre_sample)**2. + (dist_mono_sample/src_h)**2. + (dist_mono_sample/mono_h)**2. + \
    	Av_t0**2. - 2.*Av_t0*Av_t1 + Av_t1**2. ) 	# typo/missing in paper?

    Av[0,1] = Av[1,0] = 0.5*sig2fwhm**2. / ki**2. * \
        ( dist_src_mono*dist_mono_sample/src_h**2. - Av_t0*Av_t0 + Av_t0*Av_t1 )

    Av[1,1] = 0.5*sig2fwhm**2. / ki**2. * \
        ( (1./(coll_v_pre_mono))**2. + (dist_src_mono/src_h)**2. + Av_t0**2. )



	# B vector: formula 27 in [eck14]
    B = np.array([0,0,0])
    B_t0 = inv_mono_curvh / (mono_mosaic*mono_mosaic*np.abs(np.sin(thetam)))

    B[0] = sig2fwhm**2. * pos_y / ki * np.tan(thetam) * \
        ( 2.*dist_src_mono / src_w**2. + B_t0 )

    B[1] = sig2fwhm**2. * pos_y / ki * \
    ( - dist_mono_sample / (mono_w*np.abs(np.sin(thetam)))**2. + \
        B_t0 - B_t0 * inv_mono_curvh*dist_mono_sample / np.abs(np.sin(thetam)) + \
        (dist_src_mono-dist_mono_sample) / src_w**2. )



	# Bv vector: formula 39 in [eck14]
    Bv = np.array([0,0])

    Bv_t0 = inv_mono_curvv/mono_mosaic_v**2

    # typo in paper?
    Bv[0] = (-1.) *  sig2fwhm**2. * pos_z / ki * \
        ( dist_mono_sample / mono_h**2. + dist_mono_sample / src_h**2. + \
            Bv_t0 * inv_mono_curvv*dist_mono_sample - 0.5*Bv_t0 / np.abs(np.sin(thetam)) )

    # typo in paper?
    Bv[1] = (-1.) * sig2fwhm**2. * pos_z / ki * \
        ( dist_src_mono / (src_h*src_h) + 0.5*Bv_t0/np.abs(np.sin(thetam)) )




	# C scalar: formula 28 in [eck14]
    C = 0.5*sig2fwhm**2. * pos_y**2. * \
	    ( 1./src_w**2. + (1./(mono_w*np.abs(np.sin(thetam))))**2. + \
		    (inv_mono_curvh/(mono_mosaic * np.abs(np.sin(thetam))))**2. )

	# Cv scalar: formula 40 in [eck14] 
    Cv = 0.5*sig2fwhm**2. * pos_z**2. * \
        ( 1./src_h**2. + 1./mono_h**2. + (inv_mono_curvv/mono_mosaic_v)**2. )



	# z components, [eck14], equ. 42
    A[2,2] = Av[0,0] - Av[0,1]*Av[0,1]/Av[1,1]
    B[2] = Bv[0] - Bv[1]*Av(0,1)/Av(1,1)
    D = Cv - 0.25*Bv[1]/Av[1,1]


	# [eck14], equ. 54
    therefl = refl * np.sqrt(pi / Av[1,1]) # typo in paper?
    
    return [ A, B, C, D, therefl ]




#
# Eckold algorithm combining the mono and ana resolutions
#
def calc_eck(param):
    twotheta = param["twotheta"] * param["dsample_sense"]
    thetaa = param["thetaa"] * param["dana_sense"]
    thetam = param["thetam"] * param["dmono_sense"]
    ki_Q = param["angle_ki_Q"] * param["dsample_sense"]
    kf_Q = param["angle_kf_Q"] * param["dsample_sense"]


    # --------------------------------------------------------------------
    # mono/ana focus
    mono_curvh = param["mono_curvh"]
    mono_curvv = param["mono_curvv"]
    ana_curvh = param["ana_curvh"]
    ana_curvv = param["ana_curvv"]

    #if param.bMonoIsOptimallyCurvedH:
    #    mono_curvh = tl::foc_curv(param.dist_src_mono, param.dist_mono_sample, units::abs(t_real(2)*thetam), false);
    #if param.bMonoIsOptimallyCurvedV: 
    #    mono_curvv = tl::foc_curv(param.dist_src_mono, param.dist_mono_sample, units::abs(t_real(2)*thetam), true);
    #if param.bAnaIsOptimallyCurvedH: 
    #    ana_curvh = tl::foc_curv(param.dist_sample_ana, param.dist_ana_det, units::abs(t_real(2)*thetaa), false);
    #if param.bAnaIsOptimallyCurvedV: 
    #    ana_curvv = tl::foc_curv(param.dist_sample_ana, param.dist_ana_det, units::abs(t_real(2)*thetaa), true);

    inv_mono_curvh = 0.
    inv_mono_curvv = 0.
    inv_ana_curvh = 0.
    inv_ana_curvv = 0.

    if param["bMonoIsCurvedH"]:
        inv_mono_curvh = 1./mono_curvh
    if param["bMonoIsCurvedV"]: 
        inv_mono_curvv = 1./mono_curvv
    if param["bAnaIsCurvedH"]: 
        inv_ana_curvh = 1./ana_curvh
    if param["bAnaIsCurvedV"]: 
        inv_ana_curvv = 1./ana_curvv
    # --------------------------------------------------------------------


    lam = k2lam(param["ki"])

    coll_h_pre_mono = param["coll_h_pre_mono"]
    coll_v_pre_mono = param["coll_v_pre_mono"]

    if param.bGuide:
        coll_h_pre_mono = lam*param["guide_div_h"]
        coll_v_pre_mono = lam*param["guide_div_v"]



    # results
    res = {}

    res["Q_avg"] = np.array([ param["Q"], 0., 0., param["E"] ])

    # -------------------------------------------------------------------------

    # - if the instruments works in kf=const mode and the scans are counted for
    #   or normalised to monitor counts no ki^3 or kf^3 factor is needed.
    # - if the instrument works in ki=const mode the kf^3 factor is needed.

    # TODO
    #tupScFact = get_scatter_factors(param.flags, param.thetam, param.ki, param.thetaa, param.kf);
    tupScFact = [1., 1., 1.]

    dmono_refl = param["dmono_refl"] * tupScFact[0]
    dana_effic = param["dana_effic"] * tupScFact[1]
    dxsec = tupScFact[2]
    #if param.mono_refl_curve:
    #    dmono_refl *= (*param.mono_refl_curve)(param.ki)
    #if param.ana_effic_curve:
    #    dana_effic *= (*param.ana_effic_curve)(param.kf)


    #--------------------------------------------------------------------------
    # mono part

    [A, B, C, D, dReflM] = get_mono_vals(
        param["src_w"], param["src_h"],
        param["mono_w"], param["mono_h"],
        param["dist_src_mono"], param["dist_mono_sample"],
        param["ki, thetam"],
        coll_h_pre_mono, param["coll_h_pre_sample"],
        coll_v_pre_mono, param["coll_v_pre_sample"],
        param["mono_mosaic"], param["mono_mosaic_v"],
        inv_mono_curvh, inv_mono_curvv,
        param["pos_x"] , param["pos_y"], param["pos_z"],
        dmono_refl)
    #--------------------------------------------------------------------------


    #--------------------------------------------------------------------------
    # ana part
    # equ 43 in [eck14]
    pos_y2 = - param["pos_x"]*np.sin(twotheta) + param["pos_y"]*np.cos(twotheta)

    [E, F, G, H, dReflA] = get_mono_vals(
        param["det_w"], param["det_h"],
        param["ana_w"], param["ana_h"],
        param["dist_ana_det"], param["dist_sample_ana"],
        param["kf"], -thetaa,
        param["coll_h_post_ana"], param["coll_h_post_sample"],
        param["coll_v_post_ana"], param["coll_v_post_sample"],
        param["ana_mosaic"], param["ana_mosaic_v"],
        inv_ana_curvh, inv_ana_curvv,
        param["pos_x"], pos_y2, param["pos_z"],
        dana_effic)
    #--------------------------------------------------------------------------



    # equ 4 & equ 53 in [eck14]
    dE = (param["ki"]**2. - param["kf"]**2.) / (2.*param["Q"]**2.)
    kipara = param["Q"]*(0.5+dE)
    kfpara = param["Q"]-kipara
    kperp = np.sqrt(np.abs(kipara**2. - param["ki"]**2.))
    kperp *= param["dsample_sense"]



    # trafo, equ 52 in [eck14]
    T = np.identity(6)
    T[0,3] = T[1,4] = T[2,5] = -1.
    T[3,0] = 2.*ksq2E * kipara
    T[3,3] = 2.*ksq2E * kfpara
    T[3,1] = 2.*ksq2E * kperp
    T[3,4] = -2.*ksq2E * kperp
    T[4,1] = T[5,2] = (0.5 - dE)
    T[4,4] = T[5,5] = (0.5 + dE)

    Tinv = la.inv(T)


    # equ 54 in [eck14]
    Dalph_i = rotation_matrix_3d_z(-ki_Q)
    Dalph_f = rotation_matrix_3d_z(-kf_Q)
    Arot = np.dot(np.dot(np.transpose(Dalph_i), A), Dalph_i)
    Erot = np.dot(np.dot(np.transpose(Dalph_f), E), Dalph_f)

    matAE = np.zeros((6,6))
    matAE[0:3, 0:3] = Arot
    matAE[3:6, 3:6] = Erot

    # U1 matrix
    U1 = np.dot(np.dot(np.transpose(Tinv), matAE), Tinv)    # typo in paper in quadric trafo in equ 54 (top)?

    # V1 vector
    vecBF = np.zeros(6)
    vecBrot = np.dot(np.transpose(Dalph_i), B)
    vecFrot = np.dot(np.transpose(Dalph_f), F)
    vecBF[0:3] = vecBrot
    vecBF[3:6] = vecFrot
    V1 = np.dot(vecBF, Tinv)



    # --------------------------------------------------------------------------
    # integrate last 2 vars -> equs 57 & 58 in [eck14]

    U2 = quadric_proj(U1, 5);
    U = quadric_proj(U2, 4);

    V2 = quadric_proj_vec(V1, U1, 5);
    V = quadric_proj_vec(V2, U2, 4);

    W = (C + D + G + H) - 0.25*V1[5]/U1[5,5] - 0.25*V2[4]/U2[4,4]

    Z = dReflM*dReflA * np.sqrt(np.pi/np.abs(U1[5,5])) * np.sqrt(np.pi/np.abs(U2[4,4]))
    # --------------------------------------------------------------------------



    # quadratic part of quadric (matrix U)
    # careful: factor -0.5*... missing in U matrix compared to normal gaussian!
    res["reso"] = 2. * U;
    # linear (vector V) and constant (scalar W) part of quadric
    res["reso_v"] = V;
    res["reso_s"] = W;


    return res


