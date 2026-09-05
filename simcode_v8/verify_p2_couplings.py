"""NJP-120031 v8 -- P2 (Gate G2): coupling chain anchored at z_ref = 5 cm under the 1/e^2 convention.
Extends master_verify.py (review baseline). Every number entering the P2 tex edits is produced here.
Also: per-plane Lambda/Psc/SNR/kick/impulse/budget tables; CONFLICT-1 honest budget recompute;
M7 extraction from stage4_results.pkl (v7-physics values, regenerate at P7). ASCII."""
import numpy as np, pickle

hbar=1.054571817e-34; h=6.62607015e-34; c=2.99792458e8; g=9.80665
u=1.66053906660e-27; M=86.909180527*u
Gam=2*np.pi*6.07e6; Isat=25.0
d=2.5e-6; w0=2e-6; L=0.60; v0=1.5
weff=w0/np.sqrt(2)
lamD2=780.241e-9; lamD1=794.978e-9
fD2=c/lamD2; fD1=c/lamD1
zg=np.linspace(0.05,0.57,20)
vz=lambda z: np.sqrt(v0**2+2*g*z)
tz=lambda z: (vz(z)-v0)/g
tau=lambda z: w0/vz(z)
tau_ref=tau(zg[0]); tau_slit=w0/v0
k_ref=tau_ref/tau_slit
T=tz(L)

n=[0,0]
def ok(name,val,ref,tol=0.005):
    good=abs(val-ref)<=tol*abs(ref); n[0]+=good; n[1]+=1
    print(f"{'PASS' if good else 'FAIL':4s} {name:66s} calc={val:12.6g}  ms={ref:.6g}")
def hdr(s): print("="*108+f"\n{s}\n"+"="*108)

hdr("ANCHOR AND UNCHANGED OPTICS ROWS")
ok("tau_a(z_ref) (us)",tau_ref*1e6,1.11,0.005)
ok("anchor factor k = tau_ref/tau_slit",k_ref,0.8345,0.001)
D2m=(c/776.2e-9-fD2)*2*np.pi; D1m=(c/776.2e-9-fD1)*2*np.pi
D2c=(c/803.5e-9-fD2)*2*np.pi; D1c=(c/803.5e-9-fD1)*2*np.pi
ok("Delta2 meas (+THz)",D2m/2/np.pi/1e12,2.00,0.01); ok("Delta1 meas (+THz)",D1m/2/np.pi/1e12,9.12,0.005)
ok("Delta2 comp (-THz)",D2c/2/np.pi/1e12,-11.12,0.005); ok("Delta1 comp (-THz)",D1c/2/np.pi/1e12,-4.00,0.01)
bal=(2/abs(D2m)+1/abs(D1m))/(2/abs(D2c)+1/abs(D1c)); ok("balance ratio P_c/P_m (target 2.585)",bal,2.585,0.01)
I0m=2*323e-6/(np.pi*w0**2); I0c=2*835e-6/(np.pi*w0**2)
ok("I0 meas (W/m^2)",I0m,5.2e7,0.02); ok("I0 comp (W/m^2)",I0c,1.3e8,0.03)
Om_m=Gam*np.sqrt(I0m/(2*Isat)); Om_c=Gam*np.sqrt(I0c/(2*Isat))
ok("Omega0/2pi meas (GHz)",Om_m/2/np.pi/1e9,6.15,0.005); ok("Omega0/2pi comp (GHz)",Om_c/2/np.pi/1e9,9.89,0.005)
ok("zR meas (um)",np.pi*w0**2/776.2e-9*1e6,16.2,0.005); ok("zR comp (um)",np.pi*w0**2/803.5e-9*1e6,15.6,0.005)
ok("w_eff = w0/sqrt2 (um)",weff*1e6,1.41,0.005)

hdr("CONFIG E ANCHORED AT z_ref = 5 cm  (Table-2 coupling rows)")
ga_m=(Gam**2*I0m/(8*D2m*Isat))*tau_ref*100*(1+D2m/(2*D1m))
ga_c=(Gam**2*I0c/(8*D2c*Isat))*tau_ref*100*(1+D2c/(2*D1c))
ok("g_a meas (rad)",ga_m,3673,0.002); ok("g_a comp (rad)",ga_c,-3679,0.002)
om_m=2*np.pi*c/776.2e-9; om_c=2*np.pi*c/803.5e-9
Ng_m=323e-6*tau_ref/(hbar*om_m); Ng_c=835e-6*tau_ref/(hbar*om_c)
ok("N_gamma meas",Ng_m,1.40e9,0.005); ok("N_gamma comp",Ng_c,3.76e9,0.005)
phi_m=ga_m/Ng_m; phi_c=ga_c/Ng_c
ok("phi1 meas (plane-independent)",phi_m,2.62e-6,0.005); ok("phi1 comp",phi_c,-9.79e-7,0.005)
ok("phi1 cancellation: ga(z)/Ng(z) at z=0.57 == at z_ref",
   (ga_m*tau(0.57)/tau_ref)/(Ng_m*tau(0.57)/tau_ref),phi_m,1e-12)
U0=abs(ga_m)/tau_ref
ok("U0 = g_a/tau_a (rad/s, plane-independent)",U0,3.3e9,0.01)
ok("U0 plane-independence: slit-anchored value",(abs(ga_m)/k_ref)/tau_slit,U0,1e-9)
shape=lambda ww:(d/ww)**2*np.exp(-d**2/(2*ww**2))
D1f=lambda p,ww:(p**2)*shape(ww)
D1m_=D1f(phi_m,weff); D1c_=D1f(phi_c,weff)
ok("D1 meas (w_eff)",D1m_,4.48e-12,0.005); ok("D1 comp (w_eff)",D1c_,6.28e-13,0.005)
Lm=Ng_m*D1m_; Lc=Ng_c*D1c_; Leff=Lm+Lc
ok("Lambda_m",Lm,6.29e-3,0.005); ok("Lambda_c",Lc,2.36e-3,0.005); ok("Lambda_eff",Leff,8.65e-3,0.005)
V=np.exp(-Leff); D=np.sqrt(1-np.exp(-2*Leff))
ok("V",V,0.991,0.001); ok("D",D,0.131,0.005); ok("V^2+D^2",V**2+D**2,1.0,1e-9)
ok("SNR sqrt(Na Leff) @1e4 (heuristic lower bound)",np.sqrt(1e4*Leff),9.30,0.005)
ok("Na(SNR=10)",100/Leff,1.16e4,0.005)
def psc(phi,DD2,DD1):
    w2=(2/abs(DD2))/(2/abs(DD2)+1/abs(DD1)); w1=1-w2
    return abs(phi)*(w2*Gam/abs(DD2)+w1*Gam/abs(DD1))
Pm=Ng_m*psc(phi_m,D2m,D1m); Pc=Ng_c*psc(phi_c,D2c,D1c)
ok("Psc meas",Pm,0.0103,0.02); ok("Psc comp",Pc,0.0041,0.02); ok("Psc tot",Pm+Pc,0.0144,0.02)
kick_ref=np.sqrt(2/np.e)*hbar*abs(ga_m)/(M*weff)
ok("per-arm kick at z_ref (m/s), w_eff gradient",kick_ref,1.63,0.005)
sm=np.sqrt(Ng_m)/(np.sqrt(2)*abs(ga_m)); sc=np.sqrt(Ng_c)/(np.sqrt(2)*abs(ga_c))
SIG=1/np.sqrt(sm**-2+sc**-2)
ok("per-shot record noise (units hbar U0/w_eff)",SIG,6.15,0.005)
ok("per-shot SNR at <A>=1",1/SIG,0.16,0.02)
ok("m5 heuristic-vs-direct factor sqrt(2/shape)",np.sqrt(2/shape(weff)),1.75,0.005)
ok("Readout: photon flux both arms (s^-1) ~5e15",323e-6/(hbar*om_m)+835e-6/(hbar*om_c),4.64e15,0.005)
ok("Readout: N_gamma both arms per transit at z_ref",Ng_m+Ng_c,5.16e9,0.005)
ok("Psc,tot / Lambda_eff hierarchy factor",(Pm+Pc)/Leff,1.66,0.01)

hdr("FREE-SPACE REFERENCE (Eq. 17 construction, fixed Psc = 0.017) + m4 CEILING")
Del_D=2*np.pi*200e9; om_D=2*np.pi*c/780e-9
phi1_D=hbar*om_D*Gam**2/(4*np.pi*w0**2*Isat*Del_D)
ok("phi1(D) Eq.(3), w = w0 unchanged",phi1_D,2.34e-7,0.005)
ratio=phi1_D**2/(phi1_D*Gam/Del_D)
ok("phi1^2/psc1 (detuning-independent)",ratio,7.7e-3,0.005)
L1=shape(weff)*ratio*0.017
ok("Lambda_1 at Psc = 0.017 (w_eff shape)",L1,8.6e-5,0.005)
sig_res=3*(780e-9)**2/(2*np.pi)
ceil=shape(weff)*(sig_res/(2*w0**2))*0.017
ok("m4 corrected ideal ceiling sigma_res/(2 w^2)",sig_res/(2*w0**2),0.036,0.01)
ok("Lambda_1 ideal-mode-matched ceiling",ceil,4.0e-4,0.011)
ok("ceiling/actual ratio (consistency with Sec.-3 x4.7)",(sig_res/(2*w0**2))/ratio,4.7,0.005)
ok("free-space Na(SNR=10)",100/L1,1.2e6,0.04)
ok("free-space SNR @1e4",np.sqrt(1e4*L1),0.93,0.005)
taug=np.pi*w0**2/780e-9/c
gslice=(Gam**2*(2*462e-6/(np.pi*w0**2))/(2*Isat))*taug/(4*Del_D)
ok("g_slice (Config-D, unchanged)",gslice,2.28e-5,0.005)
ok("n_inst (Config-D, unchanged)",462e-6*taug/(hbar*om_D),97,0.01)
ok("first-light: g_a/30 (rad)",abs(ga_m)/30,122,0.005)
ok("first-light: Lambda_eff/30",Leff/30,2.9e-4,0.01)
ok("first-light: Na(SNR=10)",100/(Leff/30),3.5e5,0.01)

hdr("PER-PLANE TABLE (tab:couplings_planes) + IMPULSE (corrected, plane-consistent) + BUDGET")
print("    z[cm]  tau_a[us]  Leff(z)    Psc(z)   SNR@1e4  kick[m/s]  disp[um]  /w0    /weff   budget(z)")
rows={}
for kk in [0,1,2,9,19]:
    z=zg[kk]; ta=tau(z); f=ta/tau_ref
    Le=Leff*f; Ps=(Pm+Pc)*f; kz=kick_ref*f; disp=kz*ta; bud=0.1*(d/tz(z))/kz
    rows[kk]=(z,ta,Le,Ps,np.sqrt(1e4*Le),kz,disp,bud)
    print(f"    {z*100:5.2f}  {ta*1e6:8.4f}  {Le:.3e}  {Ps:.5f}  {np.sqrt(1e4*Le):6.2f}   {kz:6.3f}   {disp*1e6:7.3f}  {disp/w0:5.3f}  {disp/weff:5.3f}  {bud:.3e}")
ok("Lambda_eff at last plane",rows[19][2],4.24e-3,0.005)
ok("Psc at last plane",rows[19][3],7.05e-3,0.005)
ok("SNR@1e4 last plane",rows[19][4],6.51,0.005)
ok("kick at last plane (m/s)",rows[19][5],0.80,0.005)
ok("impulse worst plane: disp/w0",rows[0][6]/w0,0.91,0.01)
ok("impulse worst plane: disp/w_eff",rows[0][6]/weff,1.28,0.01)
ok("impulse last plane: disp/w0",rows[19][6]/w0,0.218,0.01)
budgets=[0.1*(d/tz(z))/(kick_ref*tau(z)/tau_ref) for z in zg]
ok("CONFLICT-1: honest budget = min over planes (at z = 57 cm)",min(budgets),1.42e-6,0.005)
kick_L=kick_ref*(w0/vz(L))/tau_ref
ok("CONFLICT-1: conservative z->L extrapolation",0.1*(d/T)/kick_L,1.40e-6,0.005)
print(f"NOTE plan-letter D3 value 3.2e-7 was the fixed-kinematics sqrt2 transform; honest per-plane")
print(f"     recompute (plan P2 mandate) is x{min(budgets)/3.2e-7:.1f} looser: far-plane g_a falls with tau_a while signal d/t rises.")
print(f"NOTE balanced-configuration impulse margins at eta_bal = 1e-3: barrier {hbar*U0/(0.5*M*vz(zg[0])**2)*1e-3:.2e},")
print(f"     displacement {rows[0][6]/weff*1e-3:.2e} (w_eff-relative, worst plane).")
ok("App-A barrier ratio hbar U0/(M vz(zref)^2/2) (worst plane)",hbar*U0/(0.5*M*vz(zg[0])**2),1.49,0.01)

hdr("M7 EXTRACTION FROM stage4_results.pkl  (v7-physics; regenerate at P7)")
r=pickle.load(open('simcode/sim/stage4_results.pkl','rb'))
fr7=735e-6
b={k:v/fr7*100 for k,v in r['bias_gls'].items()}
for kk,vv in b.items(): print(f"    eta,eps = {kk}: bias = {vv:.2f}% of fringe")
lo=np.mean([b[('1e-3','1e-3')],b[('1e-2','1e-4')]]); hi=b[('1e-2','1e-3')]
alpha=np.log(hi/lo)/np.log(10.0)
cross=1e-6*(5.0/lo)**(1/alpha)
print(f"    5%-of-fringe crossing (log-log interp between 1e-6 mean {lo:.2f}% and 1e-5 {hi:.1f}%): {cross:.2e}")
ok("M7: 5% crossing ~2e-6 (review value)",cross,2e-6,0.25)
ok("M7: bias at 1e-7 (%)",b[('1e-3','1e-4')],0.34,0.02)

print("\n"+"="*108); print(f"SUMMARY: {n[0]}/{n[1]} checks PASS"); print("="*108)
