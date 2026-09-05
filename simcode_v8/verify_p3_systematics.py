"""NJP-120031 v8 -- P3 (Gate G3): systematics verification. Every number entering the new
Sec. 'Environmental and apparatus systematics' is produced here. ASCII output."""
import numpy as np

h=6.62607015e-34; hbar=1.054571817e-34; kB=1.380649e-23; c=2.99792458e8
g=9.80665; u=1.66053906660e-27; M=86.909180527*u
Gam=2*np.pi*6.07e6; muB=9.2740100783e-24
d=2.5e-6; w0=2e-6; weff=w0/np.sqrt(2); L=0.60; v0=1.5
lamD2=780.241209686e-9
Om_e=7.292115e-5; lat=41.3
zg=np.linspace(0.05,0.57,20)
vz=lambda z: np.sqrt(v0**2+2*g*z); tz=lambda z:(vz(z)-v0)/g
T=tz(L); vL=vz(L)
sig=lambda z: d/tz(z)
n=[0,0]
def ok(name,val,ref,tol=0.01):
    good=abs(val-ref)<=tol*abs(ref); n[0]+=good; n[1]+=1
    print(f"{'PASS' if good else 'FAIL':4s} {name:66s} calc={val:12.6g}  ms={ref:.6g}")
def hdr(s): print("="*108+f"\n{s}\n"+"="*108)

hdr("(A) CORIOLIS  [D2; cite Lan et al. PRL 108, 090402 (2012) -- web-verified 2026-07-04]")
dvE=2*Om_e*np.cos(np.deg2rad(lat))*L
ok("eastward dv_E = 2 Om cos(lat) L (um/s)",dvE*1e6,65.7,0.005)
ok("worst-case cos(lat)=1 (um/s)",2*Om_e*L*1e6,87.5,0.005)
ok("in-plane residual @ theta = 10 mrad (um/s)",dvE*0.01*1e6,0.66,0.01)
ok("  ... as % of first-plane wiggle",dvE*0.01/sig(zg[0])*100,0.80,0.01)
ok("  ... as % of last-plane wiggle",dvE*0.01/sig(zg[19])*100,5.8,0.01)
vrec=h/(M*lamD2)
ok("transverse-velocity channel 2 Om v_rec T (um/s)",2*Om_e*vrec*T*1e6,0.20,0.03)
Ts=np.linspace(0,T,20001); zt=v0*Ts+0.5*g*Ts**2
dy=2*Om_e*np.cos(np.deg2rad(lat))*np.trapezoid(zt,Ts)
ok("out-of-plane displacement Delta y at screen (um)",dy*1e6,6.4,0.01)
ok("pre-slit eastward velocity 2 Om cos(lat) h_drop (um/s, along y)",2*Om_e*np.cos(np.deg2rad(lat))*(v0**2/(2*g))*1e6,12.6,0.01)

hdr("(B) MAGNETIC  [m10 + m11 + hyperfine-balance channel]")
a_quad=h*287.575*2*0.01*0.01/M     # |F=1,mF=0>: 287.575 Hz/G^2 per level; B=10 mG, grad=0.1 mG/cm=0.01 G/m
ok("mF=0 quadratic-Zeeman dv over T (um/s) @ {10 mG, 0.1 mG/cm}",a_quad*T*1e6,6.0e-5,0.02)
a1=muB*0.5*1e-5/M                  # |gF mF|=1/2; 1 mG/cm = 1e-5 T/m
ok("mF=+-1 accel @ 1 mG/cm (m/s^2)",a1,3.21e-4,0.005)
ok("mF=+-1 dv over T @ 1 mG/cm (um/s)",a1*T*1e6,73.5,0.005)
ok("spec {f<=5%, grad<=0.1 mG/cm}: ensemble-mean (um/s)",a1*T*0.1*0.05*1e6,0.37,0.01)
ok("  ... as % of last-plane wiggle",a1*T*0.1*0.05/sig(zg[19])*100,3.2,0.02)
ok("impurity-atom position displacement @ spec (um)",0.5*a1*0.1*T**2*1e6,0.84,0.01)
ok("hyperfine-balance channel f x 2e-3 at f = 5% (vs eta_bal 1e-3)",0.05*2e-3,1e-4,1e-9)

hdr("(C) BACKGROUND GAS  [m9]")
ok("collision prob. over T at 1e-9 mbar (%; rate 0.1/s)",0.1*T*100,2.3,0.01)
ok("collision prob. over T at 1e-10 mbar (%)",0.01*T*100,0.23,0.01)
ok("spec <= 1e-10 mbar: budget vs Lambda_eff(z_ref) ratio",0.01*T/8.65e-3,0.26,0.02)

hdr("(D) DETECTION: LIGHT-SHEET LOCALIZATION  [m8] + RECORD-ATOM PAIRING  [CONFLICT-2]")
ls=20e-6; tauc=ls/vL
ok("crossing time, 20-um sheet at v(L) (us)",tauc*1e6,5.34,0.005)
Nsc=(Gam/2)*tauc
ok("photons scattered (saturated, Gamma/2)",Nsc,102,0.01)
NA=0.5; eta_c=(1-np.cos(np.arcsin(NA)))/2*0.8*0.8
Ndet=Nsc*eta_c
ok("detected photons (NA 0.5, QE x optics 0.64)",Ndet,4.4,0.02)
sPSF=0.21*780e-9/NA
ok("diffraction PSF sigma (um)",sPSF*1e6,0.33,0.02)
ok("centroid precision sigma_PSF/sqrt(N) (um)",sPSF/np.sqrt(Ndet)*1e6,0.157,0.02)
svx=vrec*np.sqrt(Nsc/3)
ok("recoil random-walk x-blur ~ sigma_vx tauc/sqrt(3) (um)",svx*tauc/np.sqrt(3)*1e6,0.106,0.02)
ok("signal-motion blur 5.5 mm/s x tauc (um)",5.5e-3*tauc*1e6,0.029,0.03)
tot=np.sqrt((sPSF/np.sqrt(Ndet))**2+(svx*tauc/np.sqrt(3))**2+(5.5e-3*tauc)**2)
ok("total x-localization (um)  [req. <~ 2 um]",tot*1e6,0.19,0.02)
print("     arrival-time tag resolution = tauc ~ 5 us; z-thickness does not enter x.")
print("\n  PAIRING (pulsed release; transit time inferred from arrival time; residual from")
print("  longitudinal release-velocity width delta_u):  dt_pair(z) = (vL - vz)/(g vz) x delta_u")
print("    z[cm]   coeff [s/(m/s)]   dt @ du=5mm/s [us]   ceiling 1/(2e dt) [1/s]   /day")
coef={}
for kk in [0,1,2,9,19]:
    z=zg[kk]; cf=(vL-vz(z))/(g*vz(z)); coef[kk]=cf
    dt=cf*5e-3; ceil=1/(2*np.e*dt)
    print(f"    {z*100:5.2f}   {cf:10.5f}        {dt*1e6:8.1f}            {ceil:9.0f}          {ceil*86400:.2e}")
ok("pairing coeff, first plane (s per m/s)",coef[0],0.1105,0.005)
ok("pairing coeff, last plane (ms per m/s)",coef[19]*1e3,2.21,0.01)
du_opt=np.sqrt(1/(2*np.e*coef[0])*86400*0.015/8.6e7)
ok("joint-optimum delta_u (mm/s): source(du) = ceiling(du) at window top",du_opt*1e3,5.0,0.01)
thr=1/(2*np.e*coef[0]*du_opt)*86400
ok("early-plane effective throughput top (per day)",thr,2.87e7,0.01)
ok("  trim factor vs unconstrained window top (8.6e7)",8.6e7/thr,3.0,0.01)
ok("longitudinal acceptance cost du/1.5 cm/s",du_opt/0.015,0.33,0.02)
print("     clean-pairing fraction at the ceiling operating point: e^-1 = 37%; higher at lower rates.")
print("     Late planes unconstrained (z = 57 cm ceiling ~1.4e9/day). sigma_v0 (fringe smear, Sec. 2)")
print("     is set by release-height spread g sigma_z/v0 (= 3.3 mm/s at sigma_z = 0.5 mm), NOT by delta_u:")
ok("     g sigma_z/v0 at sigma_z = 0.5 mm (mm/s)  [<= 1 cm/s spec]",g*0.5e-3/v0*1e3,3.3,0.01)

hdr("(E) MULTIPASS OPTICAL BUDGET  [M3]")
Lp=0.005
Feff=(1-(1-Lp)**100)/Lp
ok("cumulative envelope at last pass (1-Lp)^99 (%)",(1-(1-Lp)**99)*100,39,0.01)
ok("effective pass sum F_eff at Lp = 0.5%",Feff,78.8,0.005)
ok("coupling shortfall F_eff/F (%)",(1-Feff/100)*100,21.2,0.01)
ok("input-power make-up factor 100/F_eff",100/Feff,1.27,0.005)
dF=((1-(1-0.0055)**100)/0.0055-Feff)/0.0005
ok("d ln F_eff / d Lp",dF/Feff,-45.0,0.01)
ok("differential-loss drift spec for eta_bal <= 1e-3: |d(DLp)|",1e-3/abs(dF/Feff),2.2e-5,0.02)
dn=6.9e-4                      # BK7 n(776.2)-n(803.5)
ok("chromatic focal shift/f = dn/(n-1) (refractive relay)",dn/0.511,1.35e-3,0.01)
ok("  -> df at f = 10 mm (um) vs zR ~ 16 um",dn/0.511*0.01*1e6,13.5,0.01)
dwreq=5e-4
ok("waist-match requirement dw/w for eta_bal_shape ~ 2 dw/w = 1e-3",2*dwreq,1e-3,1e-9)
ok("  -> focal-coincidence requirement zR sqrt(2 dw/w) (um)",16.2e-6*np.sqrt(2*dwreq)*1e6,0.51,0.02)
ratio=(2/weff**2)/(np.sqrt(2/np.e)/weff)
ok("centering coefficient |f''(0)|/|f'|_max x w_eff",ratio*weff,2.33,0.005)
ok("relative-centering requirement for eta_bal_point = 1e-3 (nm)",1e-3/ratio/1e-9,0.61,0.02)
wedge=10/3600/180*np.pi
ok("transmissive-wedge chromatic steer @ 10 arcsec, f_eff = 10 mm (nm)",dn*wedge*0.01*1e9,0.33,0.02)
print("     Vector light shift: proportional to mF -> vanishes for the |F=1, mF=0> clock state;")
print("     residual acts only on the f <= 5% Zeeman impurities and is carried in (B).")

print("\n"+"="*108); print(f"SUMMARY: {n[0]}/{n[1]} checks PASS"); print("="*108)
