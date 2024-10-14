module gamman_module
   implicit none
   public :: gamma_n, neutral_surfaces
 contains

  real function atg(s, t, p)
    ! Adiabatic temperature gradient deg C per decibar
    ! ref: Bryden, H., 1973, Deep-Sea Res., 20, 401-408
    ! Units:
    !   pressure        p        decibars
    !   temperature     t        deg Celsius (IPTS-68)
    !   salinity        s        (IPSS-78)
    !   adiabati  !   atg      deg. C/decibar
    ! Check value: atg = 3.255976e-4 C/dbar for s = 40 (IPSS-78),
    ! t = 40 deg C, p0 = 10000 decibars

    implicit none
    real, intent(in) :: s, t, p
    real :: ds

    ds = s - 35.0 

    atg = (((-2.1687d-16 * t + 1.8676d-14) * t - 4.6206d-13) * p &
          + ((2.7759d-12 * t - 1.1351d-10) * ds + ((-5.4481d-14 * t &
          + 8.733d-12) * t - 6.7795d-10) * t + 1.8741d-8)) * p &
          + (-4.2393d-8 * t + 1.8932d-6) * ds &
          + ((6.6228d-10 * t - 6.836d-8) * t + 8.5258d-6) * t + 3.5803d-5

  end function atg


 subroutine depth_ns(s, t, p, n, s0, t0, p0, sns, tns, pns)
  !   DESCRIPTION : Find the position which the neutral surface through a
  !   specified bottle intersects a neighbouring cast
  !   
  !   PRECISION : Real
  !   
  !   INPUT : s(n) array of cast salinities
  !           t(n) array of cast in situ temperatures
  !           p(n) array of cast pressures
  !           n    length of cast
  !           s0   the bottle salinity
  !           t0   the bottle in situ temperature
  !           p0   the bottle pressure
  !   
  !   OUTPUT : sns salinity of the neutral surface
  !                  intersection with the cast
  !            tns in situ temperature of the intersection
  !            pns pressure of the intersection
  !   
  !   UNITS : salinities    psu (IPSS-78)
  !           temperatures  degrees C (IPTS-68)
  !           pressures     db
  !   
  !   
  !   AUTHOR : David Jackett
  !   
  !   CREATED : June 1993
  !   
  !   REVISION : 1.1 30/6/93

    
    implicit none
    integer, parameter :: nmax = 100
    integer, intent(in) :: n
    real, intent(in) :: s(n), t(n), p(n)
    real, intent(in) :: s0, t0, p0
    real, intent(out) :: sns, tns, pns

    real :: ds, sigl, sigu, ec0, ecz0, ecz_0, eps, ez1, ez2
    real :: pc0, pc1, pc_0, p1, p2, r
    real :: sc0, tc0, ec_0 ! TODO this was not defined in the original
    integer :: k, ncr, iter, isuccess, niter
    real, dimension(nmax) :: e
    integer, parameter :: n2 = 2

    if (n > nmax) then
      print *, '\nparameter nmax in depth_ns.f < ', n, '\n'
      stop
    end if

    ncr = 0

    do k = 1, n
      call sig_vals(s0, t0, p0, s(k), t(k), p(k), sigl, sigu)
      e(k) = sigu - sigl

      if (k > 1) then
        if (e(k-1) == 0.0) then
          ncr = ncr + 1
          sns = s(k-1)
          tns = t(k-1)
          pns = p(k-1)
        elseif (e(k) * e(k-1) < 0.0) then
          ncr = ncr + 1
          pc0 = p(k-1) - e(k-1) * (p(k) - p(k-1)) / (e(k) - e(k-1))
          iter = 0
          isuccess = 0

          do while (isuccess == 0)
            iter = iter + 1
            call stp_interp(s(k-1), t(k-1), p(k-1), n2, sc0, tc0, pc0)
            call sig_vals(s0, t0, p0, sc0, tc0, pc0, sigl, sigu)
            ec0 = sigu - sigl

            p1 = (p(k-1) + pc0) / 2
            p2 = (pc0 + p(k)) / 2
            ez1 = (e(k-1) - ec0) / (pc0 - p(k-1))
            ez2 = (ec0 - e(k)) / (p(k) - pc0)
            r = (pc0 - p1) / (p2 - p1)
            ecz_0 = ez1 + r * (ez2 - ez1)

            if (iter == 1) then
              ecz0 = ecz_0
            else
              ecz0 = -(ec0 - ec_0) / (pc0 - pc_0)
              if (ecz0 == 0) ecz0 = ecz_0
            end if

            pc1 = pc0 + ec0 / ecz0
            !strategy when the iteration jumps out of the inteval

            if (pc1 <= p(k-1) .or. pc1 >= p(k)) then
              call e_solve(s, t, p, e, n, k, s0, t0, p0, sns, tns, pns, niter)
              if (pns < p(k-1) .or. pns > p(k)) then
                stop 'ERROR 1 in depth-ns.f'
              else
                isuccess = 1
              end if
            else
            !test the accuracy of the iterate

              eps = abs(pc1 - pc0)
              if (abs(ec0) <= 5.e-5 .and. eps <= 5.e-3) then
                sns = sc0
                tns = tc0
                pns = pc0
                isuccess = 1
                niter = iter
              elseif (iter > 10) then
                call e_solve(s, t, p, e, n, k, s0, t0, p0, sns, tns, pns, niter)
                isuccess = 1
              else
                pc_0 = pc0
                ec_0 = ec0
                pc0 = pc1
                isuccess = 0
              end if
            end if
          end do
        end if
      end if

      if (k == n .and. e(k) == 0.0) then
        ncr = ncr + 1
        sns = s(k)
        tns = t(k)
        pns = p(k)
      end if
    end do
    
    ! the last bottle


    if (ncr == 0) then
      sns = -99.0
      tns = -99.0
      pns = -99.0
    elseif (ncr >= 2) then
      sns = -99.2
      tns = -99.2
      pns = -99.2
    end if

end subroutine depth_ns



SUBROUTINE DERTHE(S, T, P0, DTHEDT, DTHEDS, DTHEDP)
  ! THIS SUBROUTINE USES THE BRYDEN (1973) POLYNOMIAL
  ! FOR POTENTIAL TEMPERATURE AS A FUNCTION OF S, T, P
  ! TO OBTAIN THE PARTIAL DERIVATIVES OF THETA WITH
  ! RESPECT TO T, S, P. PRESSURE IS IN DBARS.
  IMPLICIT NONE
  REAL, INTENT(IN) :: S, T, P0
  REAL, INTENT(OUT) :: DTHEDT, DTHEDS, DTHEDP

  REAL :: DS, P, PP, PPP, TT, TTT, PART
  REAL, PARAMETER :: A0 = -0.36504D-4, A1 = -0.83198D-5, A2 = +0.54065D-7
  REAL, PARAMETER :: A3 = -0.40274D-9, B0 = -0.17439D-5, B1 = +0.29778D-7
  REAL, PARAMETER :: D0 = +0.41057D-10, C0 = -0.89309D-8, C1 = +0.31628D-9
  REAL, PARAMETER :: C2 = -0.21987D-11, E0 = +0.16056D-12, E1 = -0.50484D-14

  DS = S - 35.0
  P = P0
  PP = P * P
  PPP = PP * P
  TT = T * T
  TTT = TT * T

  PART = 1.0 + P * (A1 + 2.0 * A2 * T + 3.0 * A3 * TT + DS * B1)
  DTHEDT = PART + PP * (C1 + 2.0 * C2 * T) + PPP * E1
  DTHEDS = P * (B0 + B1 * T) + PP * D0
  PART = A0 + A1 * T + A2 * TT + A3 * TTT + DS * (B0 + B1 * T)
  DTHEDP = PART + 2.0 * P * (DS * D0 + C0 + C1 * T + C2 * TT) + 3.0 * PP * (E0 + E1 * T)

  ! CHECK = T + PART * P + PP * (DS * D0 + C0 + C1 * T + C2 * TT) + PPP * (E0 + E1 * T)
  ! PRINT *, ' THE CHECK VALUE OF THETA FROM DERTHE IS = ', CHECK

  RETURN
END SUBROUTINE DERTHE


subroutine depth_scv(s, t, p, n, s0, t0, p0, sscv, tscv, pscv, nscv)
  !*********************************************************************
  ! DESCRIPTION : Find the position which the scv surface through a
  ! specified bottle intersects a neighbouring cast
  ! PRECISION : Real
  ! INPUT :
  ! s(n) array of cast salinities
  ! t(n) array of cast in situ temperatures
  ! p(n) array of cast pressures
  ! n length of cast
  ! s0 the bottle salinity
  ! t0 the bottle in situ temperature
  ! p0 the bottle pressure
  ! OUTPUT : sscv salinities of the scv surface
  ! intersections with the cast
  ! tscv temperatures of the intersections
  ! pscv pressures of the intersections
  ! nscv number of intersections
  ! UNITS : salinities psu (IPSS-78)
  ! temperatures degrees C (IPTS-68)
  ! pressures db
  ! AUTHOR : David Jackett
  ! CREATED : February 1995
  ! REVISION : 1.1 9/2/95
  !*********************************************************************

  implicit none
  integer, parameter :: n_max = 2000, nscv_max = 50
  integer :: n, ncr, k, iter, isuccess, niter, nscv
  real,intent(in) :: s0, t0, p0
  real :: sscv_tmp, tscv_tmp, pscv_tmp, sdum, sigu, sigl, ecz0, ecz_0, ec0, ec_0, eps
  real :: pc0, pc1, pc_0, p1, p2, r, ez1, ez2
  real :: sc0, tc0 ! TODO this was not defined in the original
  real,intent(in),dimension(n) :: s, t, p
  real,dimension(n) :: e
  real,intent(out),dimension(nscv_max) :: sscv, tscv, pscv

  if (n > n_max) then
    print *, '\nparameter n_max in depth-scv.f < ', n, '\n'
    stop
  end if

  ! Find the bottle pairs containing a crossing
  ncr = 0
  nscv = 0

  do k = 1, n
    sdum = svan(s0, theta(s0, t0, p0, p(k)), p(k), sigl)
    sdum = svan(s(k), t(k), p(k), sigu)
    e(k) = sigu - sigl

    if (k > 1) then
      ! An exact crossing at the k-1 bottle
      if (e(k-1) == 0.0) then
        ncr = ncr + 1
        sscv_tmp = s(k-1)
        tscv_tmp = t(k-1)
        pscv_tmp = p(k-1)
      ! A crossing between k-1 and k bottles
      elseif (e(k) * e(k-1) < 0.0) then
        ncr = ncr + 1
        ! Some Newton-Raphson iterations to find the crossing
        pc0 = p(k-1) - e(k-1) * (p(k) - p(k-1)) / (e(k) - e(k-1))
        iter = 0
        isuccess = 0

        do while (isuccess == 0)
          iter = iter + 1
          call stp_interp(s(k-1), t(k-1), p(k-1), 2, sc0, tc0, pc0)
          sdum = svan(s0, theta(s0, t0, p0, pc0), pc0, sigl)
          sdum = svan(sc0, tc0, pc0, sigu)
          ec0 = sigu - sigl

          p1 = (p(k-1) + pc0) / 2
          ez1 = (e(k-1) - ec0) / (pc0 - p(k-1))
          p2 = (pc0 + p(k)) / 2
          ez2 = (ec0 - e(k)) / (p(k) - pc0)
          r = (pc0 - p1) / (p2 - p1)
          ecz_0 = ez1 + r * (ez2 - ez1)

          if (iter == 1) then
            ecz0 = ecz_0
          else
            ecz0 = -(ec0 - ec_0) / (pc0 - pc_0)
            if (ecz0 == 0) ecz0 = ecz_0
          end if

          pc1 = pc0 + ec0 / ecz0

          ! Strategy when the iteration jumps out of the interval
          if (pc1 <= p(k-1) .or. pc1 >= p(k)) then
            call scv_solve(s, t, p, e, n, k, s0, t0, p0, sscv_tmp, tscv_tmp, pscv_tmp, niter)
            if (pscv_tmp < p(k-1) .or. pscv_tmp > p(k)) then
              stop 'ERROR 1 in depth-scv.f'
            else
              isuccess = 1
            end if
          else
            ! Otherwise, test the accuracy of the iterate
            eps = abs(pc1 - pc0)
            if (abs(ec0) <= 5.e-5 .and. eps <= 5.e-3) then
              sscv_tmp = sc0
              tscv_tmp = tc0
              pscv_tmp = pc0
              isuccess = 1
              niter = iter
            elseif (iter > 10) then
              call scv_solve(s, t, p, e, n, k, s0, t0, p0, sscv_tmp, tscv_tmp, pscv_tmp, niter)
              isuccess = 1
            else
              pc_0 = pc0
              ec_0 = ec0
              pc0 = pc1
              isuccess = 0
            end if
          end if
        end do
      end if
    end if

    ! The last bottle
    if (k == n .and. e(k) == 0.0) then
      ncr = ncr + 1
      sscv_tmp = s(k)
      tscv_tmp = t(k)
      pscv_tmp = p(k)
    end if

    ! Store multiples
    if (ncr > nscv) then
      nscv = nscv + 1
      if (nscv > nscv_max) stop 'ERROR 2 in depth-scv.f'
      sscv(nscv) = sscv_tmp
      tscv(nscv) = tscv_tmp
      pscv(nscv) = pscv_tmp
    end if
  end do

  ! No crossings
  if (nscv == 0) then
    sscv(1) = -99.0
    tscv(1) = -99.0
    pscv(1) = -99.0
  end if

end subroutine depth_scv



REAL FUNCTION EOS8D(S, T, P0, DRV)
  ! MODIFIED RCM
  ! ******************************************************
  ! SPECIFIC VOLUME ANOMALY (STERIC ANOMALY) BASED ON 1980 EQUATION
  ! OF STATE FOR SEAWATER AND 1978 PRACTICAL SALINITY SCALE.
  ! REFERENCES
  ! MILLERO, ET AL (1980) DEEP-SEA RES.,27A,255-264
  ! MILLERO AND POISSON 1981,DEEP-SEA RES.,28A PP 625-629.
  ! BOTH ABOVE REFERENCES ARE ALSO FOUND IN UNESCO REPORT 38 (1981)
  ! UNITS:      
  ! PRESSURE        P0       DECIBARS
  ! TEMPERATURE     T        DEG CELSIUS (IPTS-68)
  ! SALINITY        S        (IPSS-78)
  ! SPEC. VOL. ANA. EOS8D    M**3/KG *1.0E-8
  ! DENSITY ANA.    SIGMA    KG/M**3
  ! DRV MATRIX FORMAT
  ! 1    2     3
  ! 1   V  ,VT  ,VTT    TEMP DERIV. S,T,P
  ! 2   V0 ,VOT ,V0TT   FOR S,T,0
  ! 3   RO ,ROT ,ROTT   FOR S,T,P  DENSITY DERIV
  ! 4   K0 ,K0T ,K0TT   FOR S,T,0 SEC BULK MOD
  ! 5   A  ,AT  ,ATT
  ! 6   B  ,BT  ,BTT    BULK MOD PRESS COEFFS
  ! 7 DRDP ,K   ,DVDP   PRESSURE DERIVATIVE
  ! 8   R0S,    ,VS      SALINITY DERIVATIVES
  ! 
  ! HECK VALUE: FOR S = 40 (IPSS-78) , T = 40 DEG C, P0= 10000 DECIBARS.
  ! DR/DP                  DR/DT                 DR/DS
  ! DRV(1,7)              DRV(2,3)             DRV(1,8)
  ! 
  ! FINITE DIFFERENCE WITH 3RD ORDER CORRECTION DONE IN DOUBLE PRECSION
  ! 
  ! 3.46969238E-3       -.43311722           .705110777
  ! 
  ! EXPLICIT DIFFERENTIATION SINGLE PRECISION FORMULATION EOS80 
  ! 
  ! 3.4696929E-3        -.4331173            .7051107
  REAL, INTENT(IN) :: S, T, P0
  REAL, INTENT(OUT) :: DRV(3, 8)
  REAL :: P, SIG, SR, R1, R2, R3, R4
  REAL :: A, B, C, D, E, A1, B1, AW, BW, K, K0, KW, K35
  REAL :: R3500, DR350, V350P, SVA, SIGMA, V0, R4S, RHOS
  REAL :: RHOT, DRDT, V0T, V0S, RHOTT, V0TT, SVAN, SAL
  REAL :: BT, BTT, AT, ATT, K0S, KS, K0T, K0TT, DK, GAM, PK
  REAL :: DR35P, DVAN, VP, KT, KTT, V, V2, VS, VT, VTT, R0TT
  REAL :: DKDP, DVDP
  REAL :: RHO1, RHO2,DBDS,DADS

  DATA R3500, R4 / 1028.1063, 4.8314D-4 /
  DATA DR350 / 28.106331 /

  P = P0 / 10.0
  R3500 = 1028.1063
  SAL = S
  SR = SQRT(ABS(S))

  ! PURE WATER DENSITY AT ATMOSPHERIC PRESSURE
  ! BIGG P.H.,(1967) BR. J. APPLIED PHYSICS 8 PP 521-537.
  R1 = ((((6.536332D-9 * T - 1.120083D-6) * T + 1.001685D-4) * T - 9.095290D-3) * T + 6.793952D-2) * T - 28.263737

  ! SEAWATER DENSITY ATM PRESS. 
  ! COEFFICIENTS INVOLVING SALINITY
  ! R2 = A   IN NOTATION OF MILLERO AND POISSON 1981
  R2 = (((5.3875D-9 * T - 8.2467D-7) * T + 7.6438D-5) * T - 4.0899D-3) * T + 8.24493D-1
  R3 = (-1.6546D-6 * T + 1.0227D-4) * T - 5.72466D-3
  ! INTERNATIONAL ONE-ATMOSPHERE EQUATION OF STATE OF SEAWATER

  SIG = (R4 * S + R3 * SR + R2) * S + R1
  ! SPECIFIC VOLUME AT ATMOSPHERIC PRESSURE

  V350P = 1.0 / R3500
  SVA = -SIG * V350P / (R3500 + SIG)
  SIGMA = SIG + DR350
  DRV(1, 3) = SIGMA
  V0 = 1.0 / (1000.0 + SIGMA)
  DRV(1, 2) = V0
  ! COMPUTE DERIV WRT SALT OF RHO
  R4S = 9.6628D-4
  RHOS = R4S * S + 1.5 * R3 * SR + R2

  ! COMPUTE DERIV WRT TEMP OF RHO
  R1 = (((3.268166D-8 * T - 4.480332D-6) * T + 3.005055D-4) * T - 1.819058D-2) * T + 6.793952D-2
  R2 = ((2.155D-8 * T - 2.47401D-6) * T + 1.52876D-4) * T - 4.0899D-3
  R3 = -3.3092D-6 * T + 1.0227D-4
  RHOT = (R3 * SR + R2) * S + R1
  DRDT = RHOT
  DRV(2, 3) = RHOT
  RHO1 = 1000.0 + SIGMA
  RHO2 = RHO1 * RHO1
  V0T = -RHOT / RHO2
  !******SPECIFIC VOL. DERIV WRT S ***********
  V0S = -RHOS / RHO2
  ! COMPUTE SECOND DERIVATIVE OF RHO
  DRV(1, 8) = RHOS
  DRV(2, 2) = V0T
  ! COMPUTE SECOND DERIVATIVE OF RHO
  R1 = ((1.3072664D-7 * T - 1.3440996D-5) * T + 6.01011D-4) * T - 1.819058D-2
  R2 = (6.465D-8 * T - 4.94802D-6) * T + 1.52876D-4
  R3 = -3.3092D-6

  RHOTT = (R3 * SR + R2) * S + R1
  DRV(3, 3) = RHOTT
  V0TT = (2.0 * RHOT * RHOT / RHO1 - RHOTT) / RHO2
  ! SCALE SPECIFIC VOL. ANAMOLY TO NORMALLY REPORTED UNITS
  DRV(3, 2) = V0TT
  SVAN = SVA * 1.0E+8
  EOS8D = SVAN

  ! ******************************************************************
  ! ******  NEW HIGH PRESSURE EQUATION OF STATE FOR SEAWATER ********
  ! ******************************************************************
  ! MILLERO, ET AL , 1980 DSR 27A, PP 255-264
  ! CONSTANT NOTATION FOLLOWS ARTICLE
  !   
  ! COMPUTE COMPRESSION terms

  E = (9.1697d-10*T+2.0816d-8)*T-9.9348d-7
  BW = (5.2787d-8*T-6.12293d-6)*T+3.47718d-5
  B = BW + E*S
  ! 
  !**DERIV B WRT SALT
    DBDS=E
  !*******************
  ! CORRECT B FOR ANAMOLY BIAS CHANGE
      DRV(1,6) = B + 5.03217d-5 
  ! DERIV OF B
  BW = 1.05574d-7*T-6.12293d-6
  E = 1.83394d-9*T +2.0816d-8
  BT = BW + E*SAL
  DRV(2,6) = BT
  ! COEFFICIENTS OF A
  ! SECOND DERIV OF B
  E = 1.83394d-9
  BW = 1.05574d-7
  BTT = BW + E*SAL
  DRV(3,6) = BTT
  D = 1.91075d-4
  C = (-1.6078d-6*T-1.0981d-5)*T+2.2838d-3
  AW = ((-5.77905d-7*T+1.16092d-4)*T+1.43713d-3)*T -0.1194975
  A = (D*SR + C)*S + AW 
  ! 
  ! CORRECT A FOR ANAMOLY BIAS CHANGE
    DRV(1,5) = A + 3.3594055
  !DERIV A WRT SALT ************
    DADS=2.866125d-4*SR+C
  !*******************************
  ! DERIV OF A
  C = -3.2156d-6*T -1.0981d-5
  AW = (-1.733715d-6*T+2.32184d-4)*T+1.43713d-3
  ! 
  AT = C*SAL + AW
  DRV(2,5) = AT
  ! SECOND DERIV OF A
  C = -3.2156d-6
  AW = -3.46743d-6*T + 2.32184d-4
  ! 
  ATT = C*SAL + AW
  DRV(3,5) = ATT
  ! COEFFICIENT K0             
  B1 = (-5.3009d-4*T+1.6483d-2)*T+7.944d-2
  A1 = ((-6.1670d-5*T+1.09987d-2)*T-0.603459)*T+54.6746 
  KW = (((-5.155288d-5*T+1.360477d-2)*T-2.327105)*T +148.4206)*T-1930.06
  K0 = (B1*SR + A1)*S + KW
  ! ADD BIAS TO OUTPUT K0 VALUE
  DRV(1,4) = K0+21582.27
  !*DERIV K0 WRT SALT ************
  K0S=1.5*B1*SR+A1
  !********************************
  ! DERIV K WRT SALT   *************
  KS=(DBDS*P+DADS)*P+K0S
  !******************************
  ! DERIV OF K0
    B1 = -1.06018d-3*T+1.6483d-2
  ! APRIL 9 1984 CORRECT A1 BIAS FROM -.603457 !!!
    A1 = (-1.8501d-4*T+2.19974d-2)*T-0.603459
    KW = ((-2.0621152d-4*T+4.081431d-2)*T-4.65421)*T+148.4206
    K0T = (B1*SR+A1)*SAL + KW
    DRV(2,4) = K0T

  ! SECOND DERIV OF K0
  B1 = -1.06018d-3
  A1 = -3.7002d-4*T + 2.19974d-2
  KW = (-6.1863456d-4*T+8.162862d-2)*T-4.65421
  K0TT = (B1*SR + A1)*SAL + KW
  DRV(3,4) = K0TT

  ! EVALUATE PRESSURE POLYNOMIAL 
  ! ***********************************************
  ! K EQUALS THE SECANT BULK MODULUS OF SEAWATER
  ! DK=K(S,T,P)-K(35,0,P)
  ! K35=K(35,0,P)
  ! ***********************************************

  DK = (B * P + A) * P + K0
  K35 = (5.03217D-5 * P + 3.359406) * P + 21582.27
  GAM = P / K35
  PK = 1.0 - GAM
  SVA = SVA * PK + (V350P + SVA) * P * DK / (K35 * (K35 + DK))
  ! SCALE SPECIFIC VOL. ANAMOLY TO NORMALLY REPORTED UNITS

  SVAN = SVA * 1.0E+8
  EOS8D = SVAN
  V350P = V350P * PK

  ! ****************************************************
  ! COMPUTE DENSITY ANAMOLY WITH RESPECT TO 1000.0 KG/M**3
  ! 1) DR350: DENSITY ANAMOLY AT 35 (IPSS-78), 0 DEG. C AND 0 DECIBARS
  ! 2) DR35P: DENSITY ANAMOLY 35 (IPSS-78), 0 DEG. C ,  PRES. VARIATION
  ! 3) DVAN : DENSITY ANAMOLY VARIATIONS INVOLVING SPECFIC VOL. ANAMOLY
  ! ********************************************************************
  ! CHECK VALUE: SIGMA = 59.82037  KG/M**3 FOR S = 40 (IPSS-78),
  ! T = 40 DEG C, P0= 10000 DECIBARS.
  ! *******************************************************
  DR35P = GAM / V350P
  DVAN = SVA / (V350P * (V350P + SVA))
  SIGMA = DR350 + DR35P - DVAN
  DRV(1, 3) = SIGMA
  K = K35 + DK
  VP = 1.0 - P / K
  KT = (BT * P + AT) * P + K0T
  KTT = (BTT * P + ATT) * P + K0TT
  V = 1.0 / (SIGMA + 1000.0D0)
  DRV(1, 1) = V
  V2 = V * V
  VS = V0S * VP + V0 * P * KS / (K * K)
  RHOS = -VS / V2
  DRV(3, 8) = VS
  DRV(1, 8) = RHOS
  VT = V0T * VP + V0 * P * KT / (K * K)
  VTT = V0TT*VP+P*(2.0*V0T*KT+KTT*V0-2.0*KT*KT*V0/K)/(K*K)
  R0TT=(2.0*VT*VT/V-VTT)/V2
  DRV(3,3)=R0TT
  DRV(2,1) = VT
  DRV(3,1) = VTT
  RHOT=-VT/V2
  DRDT=RHOT
  DRV(2,3)=RHOT
  ! PRESSURE DERIVATIVE DVDP
  ! SET A & B TO UNBIASED VALUES
  A=DRV(1,5)
  B=DRV(1,6)
  DKDP = 2.0*B*P + A
  ! CORRECT DVDP TO PER DECIBAR BY MULTIPLE *.1
  DVDP = -.1*V0*(1.0 - P*DKDP/K)/K
  DRV(1,7) = -DVDP/V2
  DRV(2,7) = K
  DRV(3,7) = DVDP
END function EOS8D


SUBROUTINE EOSALL(S, T, P0, THET, SIGTHE, ALFNEW, BETNEW, GAMNEW, SOUNDV)
  ! *********************************************************************
  ! THIS PROGRAMME WRITTEN 8 JULY 1985 BY TREVOR J. McDOUGALL
  ! EOSALL STANDS FOR "EQUATION OF STATE ALL"
  ! THIS SUBROUTINE USES THE FOLLOWING FUNCTIONS WRITEN BY BOB MILLARD,
  ! - THETA(S,T,P0,PR) ; SVAN(S,T,P0,SIGTHE) ; EOSED(S,T,P0,DRV)
  ! THE NEW EXPANSION COEFFICIENT , ALFNEW , DUE TO HEAT , AND THE 
  ! CORRESPONDING SALINE CONTRACTION COEFFICIENT ARE DEFINED IN 
  ! TERMS OF THE TWO CONSERVATIVE PROPERTIES OF SEA WATER, 
  ! NAMELY POTENTIAL TEMPERATURE (REFERED TO ANY REFERENCE LEVEL)
  ! AND SALINITY. THESE COEFFICIENTS ARE DEFINED IN GILL(1982)
  ! AND HE LABELS THEM WITH DASHES, SEE HIS SECTION 3.7.4
  ! ********************************************************************
  
  IMPLICIT NONE
  REAL, INTENT(IN) :: S, T, P0
  REAL, INTENT(OUT) :: THET, SIGTHE, ALFNEW, BETNEW, GAMNEW, SOUNDV
  REAL :: PR, EDUM, ALFOLD, BETOLD, GAMOLD, SDUM
  REAL :: DTHEDT, DTHEDS, DTHEDP
  REAL :: DRV(3, 8)

  ! THE REFERENCE PRESSURE PR IS KEPT GENERAL BUT WILL BE
  ! EQUAL TO ZERO FOR ALL PERCIEVED APPLICATIONS OF NEUTRAL
  ! SURFACES.
  PR = 0.0
  THET = THETA(S, T, P0, PR)
  EDUM = EOS8D(S, T, P0, DRV)

  ! CALCULATE THE ORDINARY EXPANSION COEFFICIENTS, ALFOLD
  ! BETOLD AND THE OLD COMPRESSIBILITY, GAMOLD.
  ALFOLD = -DRV(2, 3) / (DRV(1, 3) + 1000.0)
  BETOLD = DRV(1, 8) / (DRV(1, 3) + 1000.0)
  GAMOLD = DRV(1, 7) / (DRV(1, 3) + 1000.0)
  ! CALCULATE THE SPECIFIC VOLUME ANOMALY, SVAN, AND THE
  ! SIGMA THETA, SIGTHE, BY THE FUNCTION SVAN.
  SDUM = SVAN(S, THET, PR, SIGTHE)
  CALL DERTHE(S, T, P0, DTHEDT, DTHEDS, DTHEDP)

  ALFNEW = ALFOLD / DTHEDT
  BETNEW = BETOLD + ALFNEW * DTHEDS
  GAMNEW = GAMOLD + ALFNEW * DTHEDP
  ! CHECK VALUES OF THESE NEW 'EXPANSION COEFFICIENTS'
  ! AT S=40.0,T=40.0,THET=36.8907,P0=10000.0 DBARS,
  ! ALFNEW=4395.6E-7 ; (ALFOLD=4181.1E-7)
  ! BETNEW=6646.9E-7 ; (BETOLD=6653.1E-7)
  ! GAMNEW=31.4E-7   ; (GAMOLD=32.7E-7)

  ! NOW FOR FUN WE CALCULATE THE SPEED OF SOUND. 
  ! THE IN SITU DENSITY IS (DRV(1,3)+1000.0)
  SOUNDV = SQRT(ABS(1.0E+4 / (GAMNEW * (DRV(1, 3) + 1000.0))))
  ! CHECK VALUE OF THE SPEED OF SOUND IS 
  ! AT S=40.0,T=40.0,THET=36.8907,P0=10000.0 DBARS,
  ! IS SOUNDV=1734.8 M/S.
END SUBROUTINE EOSALL



SUBROUTINE e_solve(s, t, p, e, n, k, s0, t0, p0, sns, tns, pns, iter)
  ! DESCRIPTION :	Find the zero of the e function using a 
  ! bisection method
  ! 
  ! PRECISION :		Real
  ! 
  ! INPUT :			s(n)		array of cast salinities
  ! t(n)		array of cast in situ temperatures
  ! p(n)		array of cast pressures
  ! e(n)		array of cast e values
  ! n			length of cast
  ! k			interval (k-1,k) contains the zero
  ! s0			the bottle salinity
  ! t0			the bottle in situ temperature
  ! p0			the bottle pressure
  ! 
  ! OUTPUT :		sns			salinity of the neutral surface
  ! intersection with the cast
  ! tns			in situ temperature of the intersection
  ! pns			pressure of the intersection
  ! 
  ! 
  ! UNITS :			salinities		psu (IPSS-78)
  ! temperatures	degrees C (IPTS-68)
  ! pressures		db
  ! 
  ! 
  ! AUTHOR :		David Jackett
  ! 
  ! CREATED :		June 1993
  ! 
  ! REVISION :		1.1		30/6/93
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, k
  REAL, INTENT(IN) :: s(n), t(n), p(n), e(n), s0, t0, p0
  REAL, INTENT(OUT) :: sns, tns, pns
  INTEGER, INTENT(OUT) :: iter
  REAL :: pl, el, pu, eu, pm, sm, tm, sigl, sigu, em
  INTEGER :: isuccess
  INTEGER, PARAMETER :: n2 = 2

  pl = p(k-1)
  el = e(k-1)
  pu = p(k)
  eu = e(k)

  iter = 0
  isuccess = 0

  DO WHILE (isuccess .EQ. 0)
    iter = iter + 1
    pm = (pl + pu) / 2.0

    CALL stp_interp(s(k-1), t(k-1), p(k-1), n2, sm, tm, pm)
    CALL sig_vals(s0, t0, p0, sm, tm, pm, sigl, sigu)
    em = sigu - sigl

    IF (el * em .LT. 0.0) THEN
      pu = pm
      eu = em
    ELSEIF (em * eu .LT. 0.0) THEN
      pl = pm
      el = em
    ELSEIF (em .EQ. 0.0) THEN
      sns = sm
      tns = tm
      pns = pm
      isuccess = 1
    END IF

    IF (isuccess .EQ. 0) THEN
      IF (ABS(em) .LE. 5.0E-5 .AND. ABS(pu - pl) .LE. 5.0E-3) THEN
        sns = sm
        tns = tm
        pns = pm
        isuccess = 1
      ELSEIF (iter .LE. 20) THEN
        isuccess = 0
      ELSE
        PRINT *, 'WARNING 1 in e-solve.f'
        PRINT *, iter, '  em', ABS(em), '  dp', pl, pu, ABS(pu - pl)
        sns = -99.0
        tns = -99.0
        pns = -99.0
        isuccess = 1
      END IF
    END IF
  END DO

END SUBROUTINE e_solve



subroutine gamma_n(s, t, p, n, along, alat, gamma, dg_lo, dg_hi)
  ! DESCRIPTION:
  ! Label a cast of hydrographic data at a specified 
  ! location with neutral density
  ! 
  ! PRECISION:     	Single
  ! 
  ! INPUT:
  ! s(n)        array of cast salinities
  ! t(n)        array of cast in situ temperatures
  ! p(n)        array of cast pressures
  ! n           length of cast (n=1 for single bottle)
  ! along       longitude of cast (0-360)
  ! alat        latitude of cast (-80,64)
  ! 
  ! OUTPUT:
  ! gamma(n)    array of cast gamma values
  ! dg_lo(n)    array of gamma lower error estimates
  ! dg_hi(n)    array of gamma upper error estimates
  ! 
  ! NOTE:
  ! -99.0 denotes algorithm failed
  ! -99.1 denotes input data is outside
  ! the valid range of the present
  ! equation of state
  ! 
  ! UNITS:
  ! salinity    psu (IPSS-78)
  ! temperature degrees C (IPTS-68)
  ! pressure    db
  ! gamma       kg m-3
  ! 
  ! 
  ! AUTHOR:        	David Jackett
  ! 
  ! CREATED:       	July 1993
  ! 
  ! REVISION:      	3.1     23/1/97
  implicit none
  integer, parameter :: nx=90, ny=45, nz=33, ndx=4, ndy=4
  real, parameter :: pr0=0.0, dgamma_0=0.0005, dgw_max=0.3
  integer,intent(in) :: n
  real,intent(in),dimension(n) :: s, t, p
  real,intent(out),dimension(n) :: gamma, dg_lo, dg_hi
  real, dimension(2) :: along0, alat0
  real :: along, alat
  integer :: ioce, iocean0(2, 2), n0(2, 2)
  integer :: i_min, j_min, nij, kns, itest, ij
  integer :: ialtered, i0, j0
  integer :: k
  real :: dist2_min, dist2, dx, dy, rx, ry, wsum, thk, dgamma_1, dgamma_2_l, dgamma_2_h
  real :: gw, g1_err, g2_l_err, g2_h_err, dgamma_3, rw
  real :: wt
  real, dimension(nz, 2, 2) :: s0, t0, gamma0, a0
  real, dimension(nz) :: p0
  real, dimension(4) :: gwij, wtij
  real :: sns, tns, pns  ! TODO verify these variables
  ! external :: indx, read_nc, ocean_test, depth_ns, gamma_qdr, gamma_errors, goor


  ! Error checking for input coordinates
  if (along < 0.0) then
      along = along + 360.0
      ialtered = 1
  elseif (along == 360.0) then
      along = 0.0
      ialtered = 2
  else
      ialtered = 0
  end if

  if (along < 0.0 .or. along > 360.0 .or. alat < -90.0 .or. alat > 90.0) then
      print *, 'ERROR 1 in gamma-n.f : out of oceanographic range'
      stop
  end if

  ! Initialize output arrays
  do k = 1, n
      if (s(k) < 0.0 .or. s(k) > 42.0 .or. t(k) < -2.5 .or. t(k) > 40.0 .or. p(k) < 0.0 .or. p(k) > 10000.0) then
          gamma(k) = -99.1
          dg_lo(k) = -99.1
          dg_hi(k) = -99.1
      else
          gamma(k) = 0.0
          dg_lo(k) = 0.0
          dg_hi(k) = 0.0
      end if
  end do

  ! Read records from the netCDF data file
  call read_nc(along, alat, s0, t0, p0, gamma0, a0, n0, along0, alat0, iocean0)

  ! Find the closest cast
  dist2_min = 1.0e10
  do j0 = 1, 2
      do i0 = 1, 2
          if (n0(i0, j0) /= 0) then
              dist2 = (along0(i0) - along)**2 + (alat0(j0) - alat)**2
              if (dist2 < dist2_min) then
                  i_min = i0
                  j_min = j0
                  dist2_min = dist2
              end if
          end if
      end do
  end do
  ioce = iocean0(i_min, j_min)

  ! Label the cast
  dx = abs(mod(along, real(ndx)))
  dy = abs(mod(alat + 80.0, real(ndy)))
  rx = dx / real(ndx)
  ry = dy / real(ndy)

  do k = 1, n
      if (gamma(k) /= -99.1) then
          thk = theta(s(k), t(k), p(k), pr0)
          dgamma_1 = 0.0
          dgamma_2_l = 0.0
          dgamma_2_h = 0.0
          wsum = 0.0
          nij = 0

          ! Average the gammas over the box
          do j0 = 1, 2
              do i0 = 1, 2
                  if (n0(i0, j0) /= 0) then
                      if (j0 == 1) then
                          if (i0 == 1) then
                              wt = (1.0 - rx) * (1.0 - ry)
                          elseif (i0 == 2) then
                              wt = rx * (1.0 - ry)
                          end if
                      elseif (j0 == 2) then
                          if (i0 == 1) then
                              wt = (1.0 - rx) * ry
                          elseif (i0 == 2) then
                              wt = rx * ry
                          end if
                      end if
                      wt = wt + 1.0e-6

                      call ocean_test(along, alat, ioce, along0(i0), alat0(j0), iocean0(i0, j0), p(k), itest)
                      if (itest == 0) wt = 0.0

                      call depth_ns(s0(:, i0, j0), t0(:, i0, j0), p0, n0(i0, j0), s(k), t(k), p(k), sns, tns, pns)
                      if (pns > -99.0) then
                          call indx(p0, n0(i0, j0), pns, kns)
                          call gamma_qdr(p0(kns), gamma0(kns, i0, j0), a0(kns, i0, j0), p0(kns+1), gamma0(kns+1, i0, j0), pns, gw)

                          ! Error bars
                          call gamma_errors(s0(:, i0, j0), t0(:, i0, j0), p0, &
                                            gamma0(:, i0, j0), a0(:, i0, j0), &
                                            n0(i0, j0), along0(i0), alat0(j0), &
                                            s(k), t(k), p(k), &
                                            sns, tns, pns, kns, gw, g1_err, g2_l_err, g2_h_err)
                      elseif (pns == -99.0) then
                          call goor(s0(:, i0, j0), t0(:, i0, j0), p0,&
                                    gamma0(1, i0, j0), n0(i0, j0),&
                                    s(k), t(k), p(k), gw, g1_err, g2_l_err, g2_h_err)

                          ! Adjust weight for gamma extrapolation
                          if (gw > gamma0(n0(i0, j0), i0, j0)) then
                              rw = min(dgw_max, gw - gamma0(n0(i0, j0), i0, j0)) / dgw_max
                              wt = (1.0 - rw) * wt
                          end if
                      else
                          gw = 0.0
                          g1_err = 0.0
                          g2_l_err = 0.0
                          g2_h_err = 0.0
                       end if

                      if (gw > 0.0) then
                          gamma(k) = gamma(k) + wt * gw
                          dgamma_1 = dgamma_1 + wt * g1_err
                          dgamma_2_l = max(dgamma_2_l, g2_l_err)
                          dgamma_2_h = max(dgamma_2_h, g2_h_err)
                          wsum = wsum + wt
                          nij = nij + 1
                          wtij(nij) = wt
                          gwij(nij) = gw
                      end if
                  end if
              end do
          end do

          ! The average
          if (wsum /= 0.0) then
              gamma(k) = gamma(k) / wsum
              dgamma_1 = dgamma_1 / wsum

              ! The gamma errors
              dgamma_3 = 0.0
              do ij = 1, nij
                  dgamma_3 = dgamma_3 + wtij(ij) * abs(gwij(ij) - gamma(k))
              end do
              dgamma_3 = dgamma_3 / wsum

              dg_lo(k) = max(dgamma_0, dgamma_1, dgamma_2_l, dgamma_3)
              dg_hi(k) = max(dgamma_0, dgamma_1, dgamma_2_h, dgamma_3)
          else
              gamma(k) = -99.0
              dg_lo(k) = -99.0
              dg_hi(k) = -99.0
          end if
      end if
  end do

  ! Restore original longitude if altered
  if (ialtered == 1) then
      along = along - 360.0
  elseif (ialtered == 2) then
      along = 360.0
  end if

end subroutine gamma_n


subroutine gamma_errors(s, t, p, gamma, a, n, along, alat, &
  s0, t0, p0, sns, tns, pns, kns, &
  gamma_ns, pth_error, scv_l_error, scv_h_error)
  ! DESCRIPTION :		Find the p-theta and the scv errors associated 
  ! with the basic neutral surface calculation
  ! 
  ! PRECISION :			Real
  ! 
  ! INPUT :				s(n)		array of Levitus cast salinities
  ! t(n)		array of cast in situ temperatures
  ! p(n)		array of cast pressures
  ! gamma(n)	array of cast neutral densities
  ! a(n)		array of cast quadratic coefficients
  ! n		length of cast
  ! along		longitude of Levitus cast
  ! alat		latitude of Levitus cast
  ! s0		bottle salinity
  ! t0		bottle temperature
  ! p0		bottle pressure
  ! sns		salinity of neutral surface on cast
  ! tns		temperature of neutral surface on cast
  ! pns		pressure of neutral surface on cast
  ! kns		index of neutral surface on cast
  ! gamma_ns	gamma value of neutral surface on cast
  ! 
  ! OUTPUT :			pth_error	p-theta gamma error bar
  ! scv_l_error	lower scv gamma error bar
  ! scv_h_error	upper scv gamma error bar
  ! 
  ! UNITS :				salinity	psu (IPSS-78)
  ! temperature	degrees C (IPTS-68)
  ! pressure	db
  ! gamma		kg m-3
  ! 
  ! 
  ! AUTHOR :			David Jackett
  ! 
  ! CREATED :			March 1995
  ! 
  ! REVISION :			1.1		9/3/95
implicit none

! Parameters
integer, parameter :: nscv_max = 50

! Input variables
integer,intent(in) :: n
real,intent(in) :: s(n), t(n), p(n), gamma(n), a(n)
integer ::kns
real,intent(in) :: along, alat, s0, t0, p0, sns, tns, pns, gamma_ns

! Output variables
real,intent(out) :: pth_error, scv_l_error, scv_h_error

! Local variables
real :: sscv_m(nscv_max), tscv_m(nscv_max), pscv_m(nscv_max)
real :: pr0, Tb, gamma_limit, test_limit
real :: th0, thns, sdum, rho_ns, sig_ns
real :: sig_l, sig_h, b, dp, dth, drldp, test
real :: pscv, pscv_mid, gamma_scv
integer :: kscv, nscv

! Data initialization
pr0 = 0.0
Tb = 2.7e-8
gamma_limit = 26.845
test_limit = 0.1

! p-theta error
th0 = theta(s0, t0, p0, pr0)
thns = theta(sns, tns, pns, pr0)

sdum = svan(sns, tns, pns, sig_ns)
rho_ns = 1000.0 + sig_ns

call sig_vals(s(kns), t(kns), p(kns), s(kns+1), t(kns+1), p(kns+1), sig_l, sig_h)

b = (gamma(kns+1) - gamma(kns)) / (sig_h - sig_l)

dp = pns - p0
dth = thns - th0

pth_error = rho_ns * b * Tb * abs(dp * dth) / 6.0

! scv error
scv_l_error = 0.0
scv_h_error = 0.0

if (alat <= -60.0 .or. gamma(1) >= gamma_limit) then

drldp = (sig_h - sig_l) / (rho_ns * (p(kns+1) - p(kns)))

test = Tb * dth / drldp

! Approximation
if (abs(test) <= test_limit) then
if (dp * dth >= 0.0) then
scv_h_error = (3.0 * pth_error) / (1.0 - test)
else
scv_l_error = (3.0 * pth_error) / (1.0 - test)
end if
else
! Explicit scv solution, when necessary
call depth_scv(s, t, p, n, s0, t0, p0, sscv_m, tscv_m, pscv_m, nscv)

if (nscv == 0) then
continue
else
if (nscv == 1) then
pscv = pscv_m(1)
else
pscv_mid = pscv_m((1 + nscv) / 2)
if (p0 <= pscv_mid) then
pscv = pscv_m(1)
else
pscv = pscv_m(nscv)
end if
end if

call indx(p, n, pscv, kscv)
call gamma_qdr(p(kscv), gamma(kscv), a(kscv), &
 p(kscv+1), gamma(kscv+1), pscv, gamma_scv)

if (pscv <= pns) then
scv_l_error = gamma_ns - gamma_scv
else
scv_h_error = gamma_scv - gamma_ns
end if
end if
end if
else
continue
end if

! Check for positive gamma errors
if (pth_error < 0.0 .or. scv_l_error < 0.0 .or. scv_h_error < 0.0) then
stop 'ERROR 1 in gamma-errors: negative scv error'
end if

end subroutine gamma_errors



subroutine get_lunit(lun)
  ! DESCRIPTION: Find the first FORTRAN logical unit (>=20)
  ! which is available for writing
  ! PRECISION: Real
  ! OUTPUT: lun - available logical unit

  implicit none
  integer, intent(out) :: lun
  integer :: lun0, lun1, ifound
  logical :: lv

  data lun0 /20/, lun1 /70/

  ifound = 0
  lun = lun0

  do while (lun <= lun1 .and. ifound == 0)
    inquire(unit=lun, opened=lv)
    if (lv) then
      lun = lun + 1
    else
      ifound = 1
    end if
  end do

  if (ifound == 0) stop 'ERROR 1 in get_lunit'

end subroutine get_lunit


subroutine goor(s, t, p, gamma, n, sb, tb, pb, &
                gammab,g1_err, g2_l_err, g2_h_err)
  implicit none

  ! Input parameters
  integer, intent(in) :: n
  real, intent(in) :: s(n), t(n), p(n), gamma(n)
  real, intent(in) :: sb, tb, pb

  ! Output parameters
  real, intent(out) :: gammab, g1_err, g2_l_err, g2_h_err

  ! Local variables
  real :: delt_b, delt_t, slope, pr0, Tbp
  real :: pmid, sd, sigb, sigma, s_new, t_new, e_new, s_old, t_old, e_old
  real :: sns, tns, sigl, sigu, bmid, pns, thb, thns, sig_ns, rho_ns
  real :: b, dp, dth, g2_err
  integer :: n_sth
  real :: sdum


  ! Data initialization
  delt_b = -0.1
  delt_t = 0.1
  slope = -0.14
  pr0 = 0.0
  Tbp = 2.7e-8

  ! Determine if it's bottom data
  pmid = (p(n) + pb) / 2.0
  sd = svan(s(n), theta(s(n), t(n), p(n), pmid), pmid, sigma)
  sd = svan(sb, theta(sb, tb, pb, pmid), pmid, sigb)

  ! Bottom extension
  if (sigb > sigma) then
      n_sth = 0
      s_new = s(n)
      t_new = t(n)
      e_new = sigma - sigb

      do while (sigma < sigb)
          s_old = s_new
          t_old = t_new
          e_old = e_new
          n_sth = n_sth + 1
          s_new = s(n) + n_sth * delt_b * slope
          t_new = t(n) + n_sth * delt_b
          sd = svan(s_new, theta(s_new, t_new, p(n), pmid), pmid, sigma)
          e_new = sigma - sigb
      end do

      if (sigma == sigb) then
          sns = s_new
          tns = t_new
      else
          call goor_solve(s_old, t_old, e_old, s_new, t_new, e_new, p(n), sb, tb, pb, sigb, sns, tns)
      end if

      call sig_vals(s(n-1), t(n-1), p(n-1), s(n), t(n), p(n), sigl, sigu)
      bmid = (gamma(n) - gamma(n-1)) / (sigu - sigl)

      sd = svan(s(n), t(n), p(n), sigl)
      sd = svan(sns, tns, p(n), sigu)

      gammab = gamma(n) + bmid * (sigu - sigl)
      pns = p(n)

  else
      ! Determine if it's top data
      pmid = (p(1) + pb) / 2.0
      sd = svan(s(1), theta(s(1), t(1), p(1), pmid), pmid, sigma)
      sd = svan(sb, theta(sb, tb, pb, pmid), pmid, sigb)

      ! Top extension
      if (sigb < sigma) then
          n_sth = 0
          s_new = s(1)
          t_new = t(1)
          e_new = sigma - sigb

          do while (sigma > sigb)
              s_old = s_new
              t_old = t_new
              e_old = e_new
              n_sth = n_sth + 1
              s_new = s(1)
              t_new = t(1) + n_sth * delt_t
              sd = svan(s_new, theta(s_new, t_new, p(1), pmid), pmid, sigma)
              e_new = sigma - sigb
          end do

          if (sigma == sigb) then
              sns = s_new
              tns = t_new
          else
              call goor_solve(s_new, t_new, e_new, s_old, t_old, e_old, p(1), sb, tb, pb, sigb, sns, tns)
          end if

          call sig_vals(s(1), t(1), p(1), s(2), t(2), p(2), sigl, sigu)
          bmid = (gamma(2) - gamma(1)) / (sigu - sigl)

          sd = svan(sns, tns, p(1), sigl)
          sd = svan(s(1), t(1), p(1), sigu)

          gammab = gamma(1) - bmid * (sigu - sigl)
          pns = p(1)

      else
          stop 'ERROR 1 in gamma-out-of-range.f'
      end if
  end if

  ! Error estimate
  thb = theta(sb, tb, pb, pr0)
  thns = theta(sns, tns, pns, pr0)

  sdum = svan(sns, tns, pns, sig_ns)
  rho_ns = 1000 + sig_ns

  b = bmid
  dp = pns - pb
  dth = thns - thb

  g1_err = rho_ns * b * Tbp * abs(dp * dth) / 6
  g2_err = rho_ns * b * Tbp * dp * dth / 2

  if (g2_err <= 0.0) then
      g2_l_err = -g2_err
      g2_h_err = 0.0
  else
      g2_l_err = 0.0
      g2_h_err = g2_err
  end if


end subroutine goor


subroutine goor_solve(sl, tl, el, su, tu, eu, p, s0, t0, p0, sigb, sns, tns)
  ! DESCRIPTION: Find the intersection of a potential density surface 
  ! between two bottles using a bisection method
  ! PRECISION: Real
  ! INPUT: sl, su - bottle salinities
  !        tl, tu - bottle in situ temperatures
  !        el, eu - bottle e values
  !        p - bottle pressures (the same)
  !        s0 - emanating bottle salinity
  !        t0 - emanating bottle in situ temperature
  !        p0 - emanating bottle pressure
  ! OUTPUT: sns - salinity of the neutral surface intersection with the bottle pair
  !         tns - in situ temperature of the intersection
  ! UNITS: salinities - psu (IPSS-78)
  !        temperatures - degrees C (IPTS-68)
  !        pressures - db
  ! AUTHOR: David Jackett
  ! CREATED: June 1993
  ! REVISION: 1.1 - 30/6/93

  implicit none
  real :: sl, tl, el, su, tu, eu, p, s0, t0, p0, sigb, sns, tns
  real :: rl, ru, pmid, thl, thu, rm, sm, thm, tm, sd, em
  integer :: iter, isuccess
  real :: sigma ! TODO verify this variable

  rl = 0.0
  ru = 1.0

  pmid = (p + p0) / 2.0

  thl = theta(sl, tl, p, pmid)
  thu = theta(su, tu, p, pmid)

  iter = 0
  isuccess = 0

  do while (isuccess == 0)
    iter = iter + 1

    rm = (rl + ru) / 2.0

    sm = sl + rm * (su - sl)
    thm = thl + rm * (thu - thl)

    tm = theta(sm, thm, pmid, p)

    sd = svan(sm, thm, pmid, sigma)

    em = sigma - sigb

    if (el * em < 0.0) then
      ru = rm
      eu = em
    elseif (em * eu < 0.0) then
      rl = rm
      el = em
    elseif (em == 0.0) then
      sns = sm
      tns = tm
      isuccess = 1
    end if

    if (isuccess == 0) then
      if (abs(em) <= 5.0e-5 .and. abs(ru - rl) <= 5.0e-3) then
        sns = sm
        tns = tm
        isuccess = 1
      elseif (iter <= 20) then
        isuccess = 0
      else
        print *, 'WARNING 1 in goor-solve.f'
        sns = sm
        tns = tm
        isuccess = 1
      end if
    end if
  end do

  return
end subroutine goor_solve

subroutine gamma_qdr(pl, gl, a, pu, gu, p, gamma)
  ! DESCRIPTION : Evaluate the quadratic gamma profile at a pressure
  ! between two bottles
  ! 
  ! PRECISION : Real
  ! 
  ! INPUT : pl, pu - bottle pressures
  !         gl, gu - bottle gamma values
  !         a - quadratic coefficient
  !         p - pressure for gamma value
  ! 
  ! OUTPUT : gamma - gamma value at p
  ! 
  ! UNITS : pressure - db
  !         gamma - kg m-3
  !         a - kg m-3
  ! 
  ! AUTHOR : David Jackett
  ! 
  ! CREATED : June 1993
  ! 
  ! REVISION : 1.1 30/6/93

  implicit none

  ! Input variables
  real :: pl, gl, a, pu, gu, p

  ! Output variable
  real :: gamma

  ! Local variables
  real :: p1, p2

  p1 = (p - pu) / (pu - pl)
  p2 = (p - pl) / (pu - pl)

  gamma = (a * p1 + (gu - gl)) * p2 + gl

  return
end subroutine gamma_qdr

subroutine stp_interp(s, t, p, n, s0, t0, p0)
  !
  !
  !
  !DESCRIPTION:	Interpolate salinity and in situ temperature
  !on a cast by linearly interpolating salinity
  !and potential temperature
  !
  !PRECISION:		Real
  !
  !INPUT:			s(n)	array of cast salinities
  !t(n)	array of cast in situ temperatures
  !p(n)	array of cast pressures
  !n		length of cast
  !p0		pressure for which salinity and
  !in situ temperature are required
  !
  !OUTPUT:			s0			interpolated value of salinity
  !t0			interpolated value of situ temperature
  !
  !UNITS:			salinities		psu (IPSS-78)
  !temperatures	degrees C (IPTS-68)
  !pressures		db
  !
  !
  !AUTHOR:			David Jackett
  !
  !CREATED:		June 1993
  !
  !REVISION:		1.1		30/6/93
  !
  !
  !
  implicit none

  ! Input parameters
  integer, intent(in) :: n
  real, intent(in) :: s(n), t(n), p(n), p0

  ! Output parameters
  real, intent(out) :: s0, t0

  ! Local variables
  integer :: k
  real :: r, thk, th0
  real, parameter :: pr0 = 0.0
  call indx(p, n, p0, k)
  r = (p0 - p(k)) / (p(k+1) - p(k))
  s0 = s(k) + r * (s(k+1) - s(k))
  thk = theta(s(k), t(k), p(k), pr0)
  th0 = thk + r * (theta(s(k+1), t(k+1), p(k+1), pr0) - thk)
  t0 = theta(s0, th0, pr0, p0)

end subroutine stp_interp



real function theta(s, t0, p0, pr)
  ! to compute local potential temperature at pr
  ! using bryden 1973 polynomial for adiabatic lapse rate
  ! and runge-kutta 4-th order integration algorithm.
  ! ref: bryden,h.,1973,deep-sea res.,20,401-408
  ! fofonoff,n.,1977,deep-sea res.,24,489-491
  ! units:      
  ! pressure        p0       decibars
  ! temperature     t0       deg celsius (ipts-68)
  ! salinity        s        (ipss-78)
  ! reference prs   pr       decibars
  ! potential tmp.  theta    deg celsius 
  ! checkvalue: theta= 36.89073 c,s=40 (ipss-78),t0=40 deg c,
  ! p0=10000 decibars,pr=0 decibars

  real :: s, t0, p0, pr
  real :: p, t, h, xk, q

  p = p0
  t = t0

  ! Set-up intermediate temperature and pressure variables
  h = pr - p
  xk = h * atg(s, t, p)
  t = t + 0.5 * xk
  q = xk
  p = p + 0.5 * h
  xk = h * atg(s, t, p)
  t = t + 0.29289322 * (xk - q)
  q = 0.58578644 * xk + 0.121320344 * q
  xk = h * atg(s, t, p)
  t = t + 1.707106781 * (xk - q)
  q = 3.414213562 * xk - 4.121320344 * q
  p = p + 0.5 * h
  xk = h * atg(s, t, p)
  theta = t + (xk - 2.0 * q) / 6.0

end function theta

real function svan(s, t, p0, sigma)
  ! Specific volume anomaly (steric anomaly) based on 1980 equation
  ! of state for seawater and 1978 practical salinity scale.
  ! References:
  ! Millero, et al (1980) Deep-Sea Res., 27A, 255-264
  ! Millero and Poisson 1981, Deep-Sea Res., 28A, pp 625-629.
  ! Both above references are also found in UNESCO report 38 (1981)
  ! Units:
  ! Pressure        p0       decibars
  ! Temperature     t        deg Celsius (IPTS-68)
  ! Salinity        s        (IPSS-78)
  ! Spec. vol. ana. svan     m**3/kg *1.0e-8
  ! Density ana.    sigma    kg/m**3

  real :: s, t, p0, sigma, sig
  real :: p, sr, r1, r2, r3, r4
  real :: a, b, c, d, e, a1, b1, aw, bw, k, k0, kw, k35
  real :: r3500, dr350, v350p, sva, dk, gam, pk, dr35p, dvan

  ! Data initialization
  r3500 = 1028.1063
  r4 = 4.8314e-4
  dr350 = 28.106331

  ! Convert pressure to bars and take square root of salinity
  p = p0 / 10.0
  sr = sqrt(abs(s))

  ! Pure water density at atmospheric pressure
  r1 = ((((6.536332e-9 * t - 1.120083e-6) * t + 1.001685e-4) * t - 9.095290e-3) * t + 6.793952e-2) * t - 28.263737

  ! Seawater density at atmospheric pressure
  r2 = (((5.3875e-9 * t - 8.2467e-7) * t + 7.6438e-5) * t - 4.0899e-3) * t + 8.24493e-1
  r3 = (-1.6546e-6 * t + 1.0227e-4) * t - 5.72466e-3
  sig = (r4 * s + r3 * sr + r2) * s + r1

  ! Specific volume at atmospheric pressure
  v350p = 1.0 / r3500
  sva = -sig * v350p / (r3500 + sig)
  sigma = sig + dr350
  svan = sva * 1.0e+8
  if (p == 0.0) return

  ! Compute compression terms
  e = (9.1697e-10 * t + 2.0816e-8) * t - 9.9348e-7
  bw = (5.2787e-8 * t - 6.12293e-6) * t + 3.47718e-5
  b = bw + e * s
  d = 1.91075e-4
  c = (-1.6078e-6 * t - 1.0981e-5) * t + 2.2838e-3
  aw = ((-5.77905e-7 * t + 1.16092e-4) * t + 1.43713e-3) * t - 0.1194975
  a = (d * sr + c) * s + aw
  b1 = (-5.3009e-4 * t + 1.6483e-2) * t + 7.944e-2
  a1 = ((-6.1670e-5 * t + 1.09987e-2) * t - 0.603459) * t + 54.6746
  kw = (((-5.155288e-5 * t + 1.360477e-2) * t - 2.327105) * t + 148.4206) * t - 1930.06
  k0 = (b1 * sr + a1) * s + kw

  ! Evaluate pressure polynomial
  dk = (b * p + a) * p + k0
  k35 = (5.03217e-5 * p + 3.359406) * p + 21582.27
  gam = p / k35
  pk = 1.0 - gam
  sva = sva * pk + (v350p + sva) * p * dk / (k35 * (k35 + dk))
  svan = sva * 1.0e+8
  v350p = v350p * pk

  ! Compute density anomaly
  dr35p = gam / v350p
  dvan = sva / (v350p * (v350p + sva))
  sigma = dr350 + dr35p - dvan

end function svan




subroutine indx(x, n, z, k)
  ! DESCRIPTION: Find the index of a real number in a
  ! monotonically increasing real array
  ! PRECISION: Real
  ! INPUT: x - array of increasing values
  !        n - length of array
  !        z - real number
  ! OUTPUT: k - index k if x(k) <= z < x(k+1), or
  !              n-1 if z = x(n)

  implicit none
  integer, intent(in) :: n
  real, intent(in) :: x(n), z
  integer, intent(out) :: k
  integer :: kl, ku, km

  if (x(1) < z .and. z < x(n)) then
    kl = 1
    ku = n

    do while (ku - kl > 1)
      km = (ku + kl) / 2
      if (z > x(km)) then
        kl = km
      else
        ku = km
      end if
    end do

    k = kl

    if (z == x(k + 1)) k = k + 1

  else

    if (z == x(1)) then
      k = 1
    elseif (z == x(n)) then
      k = n - 1
    else
      print *, 'ERROR 1 in indx.f : out of range'
      print *, z, n, x
    end if

  end if

end subroutine indx


subroutine neutral_surfaces(s, t, p, gamma, n, glevels, ng, sns, tns, pns, dsns, dtns, dpns)
  implicit none

  ! Parameters
  integer, parameter :: nint_max = 50
  integer :: ierr, in_error, nint, int_middle, i_int, k, ig, lun
  integer :: int(nint_max)
  integer,intent(in) :: n, ng
  real,intent(in) :: s(n), t(n), p(n), gamma(n), glevels(ng)
  real, intent(out):: sns(ng), tns(ng), pns(ng), dsns(ng), dtns(ng), dpns(ng)
  real :: alfa_l, beta_l, alfa_u, beta_u, alfa_mid, beta_mid
  real :: rhomid, thl, thu, dels, delth, pl, pu, delp, delp2, bmid
  real :: a, b, c, q, pmid, thdum, sthdum, alfa, beta, gdum, sdum
  real :: smid, tmid, sigmid, bden, pns1, pns2, rg, sd
  real :: sns_top, tns_top, pns_top, sns_middle, tns_middle, pns_middle
  real, parameter :: pr0 = 0.0, ptol = 1.0e-3
  integer, parameter :: n2 = 2

  ! Detect error condition
  in_error = 0
  do k = 1, n
      if (gamma(k) <= 0.0) in_error = 1
  end do

  if (in_error == 1) stop '\nERROR 1 in neutral-surfaces.f : missing gamma value'

  ! Loop over the surfaces
  ierr = 0
  do ig = 1, ng
      nint = 0
      do k = 1, n-1
          if (min(gamma(k), gamma(k+1)) <= glevels(ig) .and. glevels(ig) <= max(gamma(k), gamma(k+1))) then
              nint = nint + 1
              if (nint > nint_max) stop 'ERROR 2 in neutral-surfaces.f'
              int(nint) = k
          end if
      end do

      if (nint == 0) then
          sns(ig) = -99.0
          tns(ig) = -99.0
          pns(ig) = -99.0
          dsns(ig) = 0.0
          dtns(ig) = 0.0
          dpns(ig) = 0.0
      else
          if (mod(nint, 2) == 0 .and. int(1) > n/2) then
              int_middle = (nint + 2) / 2
          else
              int_middle = (nint + 1) / 2
          end if

          do i_int = 1, nint
              k = int(i_int)
              pmid = (p(k) + p(k+1)) / 2.0

              call eosall(s(k), t(k), p(k), thdum, sthdum, alfa, beta, gdum, sdum)
              alfa_l = alfa
              beta_l = beta
              call eosall(s(k+1), t(k+1), p(k+1), thdum, sthdum, alfa, beta, gdum, sdum)
              alfa_u = alfa
              beta_u = beta

              alfa_mid = (alfa_l + alfa_u) / 2.0
              beta_mid = (beta_l + beta_u) / 2.0

              call stp_interp(s(k), t(k), p(k), n2, smid, tmid, pmid)

              sd=svan(smid, tmid, pmid, sigmid)
              rhomid = 1000.0 + sigmid

              thl = theta(s(k), t(k), p(k), pr0)
              thu = theta(s(k+1), t(k+1), p(k+1), pr0)

              dels = s(k+1) - s(k)
              delth = thu - thl

              pl = p(k)
              pu = p(k+1)
              delp = pu - pl
              delp2 = delp * delp

              bden = rhomid * (beta_mid * dels - alfa_mid * delth)
              if (abs(bden) <= 1.0e-6) bden = 1.0e-6

              bmid = (gamma(k+1) - gamma(k)) / bden

              a = dels * (beta_u - beta_l) - delth * (alfa_u - alfa_l)
              a = (a * bmid * rhomid) / (2 * delp2)

              b = dels * (pu * beta_l - pl * beta_u) - delth * (pu * alfa_l - pl * alfa_u)
              b = (b * bmid * rhomid) / delp2

              c = dels * (beta_l * (pl - 2.0 * pu) + beta_u * pl) - delth * (alfa_l * (pl - 2.0 * pu) + alfa_u * pl)
              c = gamma(k) + (bmid * rhomid * pl * c) / (2 * delp2)
              c = c - glevels(ig)

              if (a /= 0.0 .and. bden /= 1.0e-6) then
                  q = -(b + sign(1.0, b) * sqrt(b * b - 4 * a * c)) / 2.0
                  pns1 = q / a
                  pns2 = c / q

                  if (pns1 >= p(k) - ptol .and. pns1 <= p(k+1) + ptol) then
                      pns(ig) = min(p(k+1), max(pns1, p(k)))
                  elseif (pns2 >= p(k) - ptol .and. pns2 <= p(k+1) + ptol) then
                      pns(ig) = min(p(k+1), max(pns2, p(k)))
                  else
                      stop 'ERROR 3 in neutral-surfaces.f'
                  end if
              else
                  rg = (glevels(ig) - gamma(k)) / (gamma(k+1) - gamma(k))
                  pns(ig) = p(k) + rg * (p(k+1) - p(k))
              end if

              call stp_interp(s, t, p, n, sns(ig), tns(ig), pns(ig))

              if (nint > 1) then
                  if (ierr == 0) then
                      ierr = 1
                      call get_lunit(lun)
                      open(lun, file='ns-multi.dat', status='unknown')
                  end if

                  if (i_int == 1) write(lun, *) ig, nint
                  write(lun, *) sns(ig), tns(ig), pns(ig)

                  if (i_int == 1) then
                      sns_top = sns(ig)
                      tns_top = tns(ig)
                      pns_top = pns(ig)
                  elseif (i_int == int_middle) then
                      sns_middle = sns(ig)
                      tns_middle = tns(ig)
                      pns_middle = pns(ig)
                  elseif (i_int == nint) then
                      if ((pns_middle - pns_top) > (pns(ig) - pns_middle)) then
                          dsns(ig) = sns_middle - sns_top
                          dtns(ig) = tns_middle - tns_top
                          dpns(ig) = pns_middle - pns_top
                      else
                          dsns(ig) = sns(ig) - sns_middle
                          dtns(ig) = tns(ig) - tns_middle
                          dpns(ig) = pns(ig) - pns_middle
                      end if
                      sns(ig) = sns_middle
                      tns(ig) = tns_middle
                      pns(ig) = pns_middle
                  end if
              else
                  dsns(ig) = 0.0
                  dtns(ig) = 0.0
                  dpns(ig) = 0.0
              end if
          end do
      end if
  end do

  if (ierr == 1) close(lun)

end subroutine neutral_surfaces


subroutine ocean_test(x1, y1, io1, x2, y2, io2, z, itest)
  ! DESCRIPTION: Test whether two locations are connected by ocean
  ! 
  ! PRECISION: Real
  ! 
  ! INPUT: x1 - longitude of first location
  !        y1 - latitude of first location
  !        io1 - ocean of first location
  !        x2 - longitude of second location
  !        y2 - latitude of second location
  !        io2 - ocean of second location
  !        z - depth of connection
  ! 
  ! OUTPUT: itest - success of connection
  ! 
  ! AUTHOR: David Jackett
  ! 
  ! CREATED: June 1994
  ! 
  ! REVISION: 1.1 7/7/94

  implicit none

  ! Input variables
  real :: x1, y1, x2, y2, z
  integer :: io1, io2

  ! Output variable
  integer :: itest

  ! Local variables
  real :: y, em1, c1, em2, c2
  integer :: isj1, isj2
  real, dimension(3) :: x_js = [129.87, 140.37, 142.83]
  real, dimension(3) :: y_js = [32.75, 37.38, 53.58]

  y = (y1 + y2) / 2.0

  ! Same ocean talks
  if (io1 == io2) then
    itest = 1
    return
  elseif (y <= -20.0) then
    ! Land of South America doesn't talk
    if (y >= -48.0 .and. (io1 * io2) == 12) then
      itest = 0
    else
      itest = 1
    end if
  elseif ((io1 == 1 .or. io1 == 2) .and. (io2 == 1 .or. io2 == 2)) then
    ! Pacific talks
    itest = 1
  elseif ((io1 == 3 .or. io1 == 4) .and. (io2 == 3 .or. io2 == 4)) then
    ! Indian talks
    itest = 1
  elseif ((io1 == 5 .or. io1 == 6) .and. (io2 == 5 .or. io2 == 6)) then
    ! Atlantic talks
    itest = 1
  elseif ((io1 * io2) == 8 .and. z <= 1200.0 .and. x1 >= 124.0 .and. x1 <= 132.0 .and. x2 >= 124.0 .and. x2 <= 132.0) then
    ! Indonesian throughflow
    itest = 1
  else
    ! Anything else doesn't
    itest = 0
  end if

  ! Exclude Japan Sea from talking
  if ((x_js(1) <= x1 .and. x1 <= x_js(3) .and. y_js(1) <= y1 .and. y1 <= y_js(3)) .or. &
      (x_js(1) <= x2 .and. x2 <= x_js(3) .and. y_js(1) <= y2 .and. y2 <= y_js(3))) then

    em1 = (y_js(2) - y_js(1)) / (x_js(2) - x_js(1))
    c1 = y_js(1) - em1 * x_js(1)

    em2 = (y_js(3) - y_js(2)) / (x_js(3) - x_js(2))
    c2 = y_js(2) - em2 * x_js(2)

    if ((y1 - em1 * x1 - c1) >= 0.0 .and. (y1 - em2 * x1 - c2) >= 0.0) then
      isj1 = 1
    else
      isj1 = 0
    end if

    if ((y2 - em1 * x2 - c1) >= 0.0 .and. (y2 - em2 * x2 - c2) >= 0.0) then
      isj2 = 1
    else
      isj2 = 0
    end if

    if (isj1 == isj2) then
      itest = 1
    else
      itest = 0
    end if
  end if

  ! Exclude Antarctic tip
  if ((io1 * io2) == 12 .and. y < -60.0) itest = 0

  return
end subroutine ocean_test



integer function ilong(alng, ndx)
  real :: alng
  integer :: ndx
  ilong = int(alng / ndx + 1)
end function ilong


integer function jlat(alt, ndy)
  real :: alt
  integer :: ndy
  jlat = int((88.0 + alt) / ndy + 1)
end function jlat


subroutine read_nc(along, alat, s0, t0, p0, gamma0, a0, n0, along0, alat0, iocean0)
  implicit none

  ! Parameters
  integer, parameter :: nx=90, ny=45, nz=33, ndx=4, ndy=4

  ! Input
  real, intent(in) :: along, alat

  ! Output
  real, intent(out) :: s0(nz,2,2), t0(nz,2,2), p0(nz), gamma0(nz,2,2), a0(nz,2,2)
  integer, intent(out) :: n0(2,2), iocean0(2,2)
  real, intent(out) :: along0(2), alat0(2)

  ! Local variables
  integer :: n(nx,ny), i0, j0, i, j, k, lun, irec, jrec, krec, lun1
  integer :: iocean(nx,ny)
  real :: along_s(nx), alat_s(ny), p0_s(nz)
  real :: s0_s(nz,2,2), t0_s(nz,2,2), gamma0_s(nz,2,2), a0_s(nz,2,2)
  real :: along_d(nx), alat_d(ny)
  real :: dx, dy

  ! Save variables
  save along_d, alat_d, p0_s, n, iocean, i0, j0

  ! Data initialization
  data i0 /1/, j0 /1/


  ! ! Calculate indices
  ! i0 = int(along / ndx + 1)
  ! j0 = int((88 + alat) / ndy + 1)


  ! Only read when necessary
  dx = along - along_d(i0)
  dy = alat - alat_d(j0)
  

  if (dx < 0.0 .or. dx >= 4.0 &
        .or. dy < 0.0 .or. dy >= 4.0 .or. (i0 == 1 .and. j0 == 1)) then
    ! Read the 'llp.fdt' file
    call get_lunit(lun)
    
    if (i0 == 1 .and. j0 == 1) then
      open(lun, file='/home/otel/Desktop/deletar/pygamman/pygamman/fortrandata/llp.fdt',&
      status='old', form='unformatted')
      read(lun) along_s, alat_s, p0_s, n, iocean
      close(lun)

      do i = 1, nx
        along_d(i) = along_s(i)
      end do

      do j = 1, ny
        alat_d(j) = alat_s(j)
      end do
    end if

    do k = 1, nz
      p0(k) = p0_s(k)
    end do

    ! Read the appropriate records from 'stga.fdt'
    i0 = int(ilong(along,ndx))
    j0 = int(jlat(alat,ndy))

    if (i0 == nx + 1) i0 = 1

    along0(1) = along_d(i0)
    alat0(1) = alat_d(j0)
    alat0(2) = alat0(1) + ndy

    ! open(lun, file='/home/guido/anaconda2/lib/python2.7/site-packages/pygamman/ocndata/stga.fdt', status='old', access='direct', recl=528, form='unformatted')
    open(lun1, file='/home/otel/Desktop/deletar/pygamman/pygamman/fortrandata/stga.fdt',&
      status='old', access='direct', recl=528, form='unformatted')  !TODO review this line

    if (i0 < nx) then
      along0(2) = along0(1) + ndx

      irec = i0
      jrec = j0
      krec = irec + (jrec - 1) * nx
      read(lun1, rec=krec) (s0_s(k,1,1), k=1,nz), (t0_s(k,1,1), k=1,nz), (gamma0_s(k,1,1), k=1,nz), (a0_s(k,1,1), k=1,nz)
      


      irec = i0 + 1
      jrec = j0
      krec = irec + (jrec - 1) * nx
      read(lun1, rec=krec) (s0_s(k,2,1), k=1,nz), (t0_s(k,2,1), k=1,nz), (gamma0_s(k,2,1), k=1,nz), (a0_s(k,2,1), k=1,nz)

      irec = i0
      jrec = j0 + 1
      krec = irec + (jrec - 1) * nx
      read(lun1, rec=krec) (s0_s(k,1,2), k=1,nz), (t0_s(k,1,2), k=1,nz), (gamma0_s(k,1,2), k=1,nz), (a0_s(k,1,2), k=1,nz)

      irec = i0 + 1
      jrec = j0 + 1
      krec = irec + (jrec - 1) * nx
      read(lun1, rec=krec) (s0_s(k,2,2), k=1,nz), (t0_s(k,2,2), k=1,nz), (gamma0_s(k,2,2), k=1,nz), (a0_s(k,2,2), k=1,nz)
      

    elseif (i0 == nx) then
      along0(2) = 0.0

      irec = i0
      jrec = j0
      krec = irec + (jrec - 1) * nx
      read(lun1, rec=krec) (s0_s(k,1,1), k=1,nz), (t0_s(k,1,1), k=1,nz), (gamma0_s(k,1,1), k=1,nz), (a0_s(k,1,1), k=1,nz)

      irec = i0
      jrec = j0 + 1
      krec = irec + (jrec - 1) * nx
      read(lun1, rec=krec) (s0_s(k,1,2), k=1,nz), (t0_s(k,1,2), k=1,nz), (gamma0_s(k,1,2), k=1,nz), (a0_s(k,1,2), k=1,nz)

      irec = 1
      jrec = j0
      krec = irec + (jrec - 1) * nx
      read(lun1, rec=krec) (s0_s(k,2,1), k=1,nz), (t0_s(k,2,1), k=1,nz), (gamma0_s(k,2,1), k=1,nz), (a0_s(k,2,1), k=1,nz)

      irec = 1
      jrec = j0 + 1
      krec = irec + (jrec - 1) * nx
      read(lun1, rec=krec) (s0_s(k,2,2), k=1,nz), (t0_s(k,2,2), k=1,nz), (gamma0_s(k,2,2), k=1,nz), (a0_s(k,2,2), k=1,nz)
    end if

    ! Get the depth and ocean information
    do j = 1, 2
      do i = 1, 2
        
        n0(i,j) = n(ilong(along0(i), ndx), jlat(alat0(j), ndy))
        iocean0(i,j) = iocean(ilong(along0(i), ndx), jlat(alat0(j), ndy))
      end do
    end do

    ! Transfer the data
    do j = 1, 2
      do i = 1, 2
        do k = 1, nz
          s0(k,i,j) = s0_s(k,i,j)
          t0(k,i,j) = t0_s(k,i,j)
          gamma0(k,i,j) = gamma0_s(k,i,j)
          a0(k,i,j) = a0_s(k,i,j)
        end do
      end do
    end do

    close(lun1)
  end if

end subroutine read_nc


SUBROUTINE scv_solve(s, t, p, e, n, k, s0, t0, p0, sscv, tscv, pscv, iter)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, k
  REAL, INTENT(IN) :: s(n), t(n), p(n), e(n)
  REAL, INTENT(IN) :: s0, t0, p0
  REAL, INTENT(OUT) :: sscv, tscv, pscv
  INTEGER, INTENT(OUT) :: iter

  REAL :: pl, el, pu, eu, pm, sm, tm, sdum, sigl, sigu, em
  INTEGER :: isuccess
  INTEGER, PARAMETER :: n2 = 2

  pl = p(k-1)
  el = e(k-1)
  pu = p(k)
  eu = e(k)

  iter = 0
  isuccess = 0

  DO WHILE (isuccess == 0)
    iter = iter + 1
    pm = (pl + pu) / 2.0

    CALL stp_interp(s(k-1), t(k-1), p(k-1), n2, sm, tm, pm)

    sdum = svan(s0, theta(s0, t0, p0, pm), pm, sigl)
    sdum = svan(sm, tm, pm, sigu)
    em = sigu - sigl

    IF (el * em < 0.0) THEN
      pu = pm
      eu = em
    ELSEIF (em * eu < 0.0) THEN
      pl = pm
      el = em
    ELSEIF (em == 0.0) THEN
      sscv = sm
      tscv = tm
      pscv = pm
      isuccess = 1
    END IF

    IF (isuccess == 0) THEN
      IF (ABS(em) <= 5.0E-5 .AND. ABS(pu - pl) <= 5.0E-3) THEN
        sscv = sm
        tscv = tm
        pscv = pm
        isuccess = 1
      ELSEIF (iter <= 20) THEN
        isuccess = 0
      ELSE
        PRINT *, 'WARNING 1 in scv-solve.f'
        PRINT *, iter, '  em', ABS(em), '  dp', pl, pu, ABS(pu - pl)
        sscv = -99.0
        tscv = -99.0
        pscv = -99.0
        isuccess = 1
      END IF
    END IF
  END DO

END SUBROUTINE scv_solve



SUBROUTINE sig_vals(s1, t1, p1, s2, t2, p2, sig1, sig2)
  IMPLICIT NONE
  REAL, INTENT(IN) :: s1, t1, p1, s2, t2, p2
  REAL, INTENT(OUT) :: sig1, sig2
  REAL :: pmid, sd

  pmid = (p1 + p2) / 2.0

  sd = svan(s1, theta(s1, t1, p1, pmid), pmid, sig1)
  sd = svan(s2, theta(s2, t2, p2, pmid), pmid, sig2)

END SUBROUTINE sig_vals



end module gamman_module
