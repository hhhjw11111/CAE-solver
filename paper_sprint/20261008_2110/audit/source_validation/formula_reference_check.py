"""Pure mathematical reference checks; NEVER a production solver validation."""
from pathlib import Path
import json,csv,math,hashlib,datetime
out=Path(__file__).absolute().parents[2]/'validation/formula_reference'
out.mkdir(parents=True,exist_ok=False)
contract={'classification':'PURE_MATHEMATICAL_REFERENCE_NOT_PRODUCTION_VALIDATION','A_values':[0,1,10,40,-40],
 'normalization':'x/L in [0,1]; n/nL; nR/nL=2; current divided by q D nL/L; A=mu E L/D',
 'endpoint_absolute_tolerance':5e-14,'constant_current_scaled_tolerance':5e-13,
 'power_identity_scaled_tolerance':2e-10,'zero_field_limit_absolute_tolerance':2e-8,
 'zero_limit_A':[-1e-8,1e-8],'simpson_panels':4096,
 'thermal_identity_scaled_tolerance':5e-14,
 'sampling':[i/100 for i in range(101)],'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(out/'frozen_contract.json').write_text(json.dumps(contract,indent=2),encoding='utf-8')
def state(A,x):
    if A==0:return 1+x,1.,1.
    if A>0:
        den=-math.expm1(-A)
        n=1-math.expm1(-A*x)/den
        grad=A*math.exp(-A*x)/den
        exact=A*(1+1/den)
    else:
        b=-A;den=-math.expm1(-b)
        ratio=math.exp(b*(x-1))*(-math.expm1(-b*x))/den
        n=1+ratio
        grad=b*math.exp(b*(x-1))/den
        exact=b*(2*math.exp(-b)-1)/den
    return n,grad,exact
def simpson(fun,N):
    return (fun(0)+fun(1)+sum((4 if i%2 else 2)*fun(i/N) for i in range(1,N)))/(3*N)
rows=[];checks=[]
for A in contract['A_values']:
    maxcur=0
    for x in contract['sampling']:
        n,g,j=state(A,x);jr=A*n+g
        rows.append({'A':A,'x_over_L':x,'n_over_nL':n,'dn_normalized':g,'J_drift_normalized':A*n,'J_diff_normalized':g,'J_total_normalized':jr,'J_exact_normalized':j})
        maxcur=max(maxcur,abs(jr-j)/max(1,abs(j)))
    endpoints=max(abs(state(A,0)[0]-1),abs(state(A,1)[0]-2))
    j=state(A,0)[2]
    pd=simpson(lambda x:A*A*state(A,x)[0],contract['simpson_panels'])
    pdiff=A
    pwork=A*j
    pterm=-((-A)*j-0*j)
    pe=abs(pd+pdiff-pwork)/max(1,abs(pwork),abs(pd))
    checks.append({'A':A,'endpoint_error':endpoints,'maximum_constant_current_scaled_error':maxcur,'P_drift_normalized':pd,'P_diffusion_normalized':pdiff,'P_work_normalized':pwork,'P_terminal_normalized':pterm,'power_identity_scaled_error':pe,'passed':endpoints<=contract['endpoint_absolute_tolerance'] and maxcur<=contract['constant_current_scaled_tolerance'] and pe<=contract['power_identity_scaled_tolerance']})
with (out/'drift_diffusion_reference.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
limits=[]
for A in contract['zero_limit_A']:
    nerr=max(abs(state(A,x)[0]-(1+x)) for x in contract['sampling'])
    jerr=abs(state(A,.5)[2]-1)
    limits.append({'A':A,'max_density_zero_limit_error':nerr,'current_zero_limit_error':jerr,'passed':max(nerr,jerr)<=contract['zero_field_limit_absolute_tolerance']})
# Actual input constants used only in a separately labelled constant-ohmic reference.
q=1.602176634e-19;n0=6e22;mu0=.0205;k=60.;L=20e-6;V=3.
sigma=q*n0*mu0;thermal=[]
for x in contract['sampling']:
    psi=V*x;theta=sigma*psi*(V-psi)/(2*k)
    thermal.append({'x_over_L':x,'psi_V':psi,'theta_reference_K':theta,'classification':'CONSTANT_OHMIC_REFERENCE_NOT_SOLVER_RESULT'})
thermal_res=abs(-k*(-sigma*V*V/(k*L*L))-sigma*(V/L)**2)/max(1,sigma*(V/L)**2)
with (out/'thermal_reference_3V.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=thermal[0]);w.writeheader();w.writerows(thermal)
counter=[]
for x in contract['sampling']:
    n=math.exp(-x);jd=n;jdiff=-n
    counter.append({'x_over_L':x,'n_over_nL':n,'J_drift':jd,'J_diffusion':jdiff,'J_total':jd+jdiff,'Q_drift':n,'J_total_dot_E':0.})
with (out/'local_equilibrium_heat_counterexample.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=counter[0]);w.writeheader();w.writerows(counter)
result={'classification':contract['classification'],'cases':checks,'zero_limit':limits,
 'thermal_reference':{'sigma_S_per_m':sigma,'k_W_per_mK':k,'L_m':L,'V_V':V,'theta_max_K':sigma*V*V/(8*k),'pde_scaled_residual':thermal_res,'scope':'homogeneous isotropic constant-ohmic model, end temperatures fixed, remaining boundaries insulating; not ZnO coupled solution'},
 'counterexample':'At A=1 and n=exp(-x), drift and diffusion currents cancel exactly, while Q_drift=exp(-x)>0; current work is zero.',
 'formula_check_passed':all(c['passed'] for c in checks+limits) and thermal_res<=contract['thermal_identity_scaled_tolerance'],
 'production_solver_validation':False,'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(out/'formula_reference_summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
