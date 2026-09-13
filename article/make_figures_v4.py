from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT=Path('/mnt/data/revised_lognormal_mixture_article_v4_package/figures')
OUT.mkdir(exist_ok=True)
FX10=Path('/mnt/data/work_fx_bal/VolatilitySmileTwoModels_BalancedBoundary_Equities_SigmaFloor10pct/fx10_runs')
FX02=Path('/mnt/data/work_fx_bal/VolatilitySmileTwoModels_BalancedBoundary_Equities_SigmaFloor10pct/fx2_runs')
EQ10=Path('/mnt/data/work_equity_eta10/VolatilitySmileTwoModels_RichEurostox30_BetterPlots/eta10_runs')
EQ02=Path('/mnt/data/work_equity_eta2/VolatilitySmileTwoModels_RichEurostox30_BetterPlots/eta2_runs')
STRESS=Path('/mnt/data/work_stress/VolatilitySmileTwoModels_EquityBoundedKernels_RightWing/output')

plt.rcParams.update({'font.size':9.3,'axes.titlesize':10.5,'axes.labelsize':9.5,'legend.fontsize':7.8,'figure.dpi':150})

def save(fig,name):
    fig.tight_layout()
    fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight')
    fig.savefig(OUT/f'{name}.png',dpi=190,bbox_inches='tight')
    plt.close(fig)

def fxdata(root,mode,file):
    return pd.read_csv(root/mode/f'{mode}_{file}.csv')

def fx_base(mode,mat,name,title):
    sm=fxdata(FX10,mode,'continuous_smiles'); de=fxdata(FX10,mode,'densities'); ma=fxdata(FX10,mode,'market_vs_models')
    s=sm[sm.maturity.astype(str)==mat].sort_values('k'); d=de[de.maturity.astype(str)==mat].sort_values('k'); q=ma[ma.maturity.astype(str)==mat].sort_values('k')
    klo=float(q.k.min()); khi=float(q.k.max()); w=khi-klo
    xmin=max(float(s.k.min()),klo-0.85*w); xmax=min(float(s.k.max()),khi+1.25*w)
    fig,axs=plt.subplots(1,2,figsize=(10.6,4.0))
    ax=axs[0]; ax.plot(s.k,100*s.model1_lnmix_iv,lw=2,label=r'Mixture ($\eta=10\%$)'); ax.plot(s.k,100*s.model2_cubic_iv,lw=1.6,ls='--',label='Cubic total-variance spline'); ax.scatter(q.k,100*q.market_iv,s=11,marker='s',zorder=5,label='Market quotes'); ax.axvline(klo,ls=':',lw=.8);ax.axvline(khi,ls=':',lw=.8);ax.set_xlim(xmin,xmax);ax.set_xlabel('Moneyness $K/F$');ax.set_ylabel('Implied volatility (%)');ax.set_title(title+' - smile');ax.grid(alpha=.2);ax.legend()
    ax=axs[1]; ax.plot(d.k,d.model1_density,lw=2,label=r'Mixture density ($\eta=10\%$)');ax.plot(d.k,d.model2_density,lw=1.6,ls='--',label='Cubic-spline density');ax.axhline(0,lw=.7);ax.axvline(klo,ls=':',lw=.8);ax.axvline(khi,ls=':',lw=.8);ax.set_xlim(xmin,xmax);ax.set_xlabel('Moneyness $K/F$');ax.set_ylabel('Normalized risk-neutral density');ax.set_title(title+' - density');ax.grid(alpha=.2);ax.legend();save(fig,name)

fx_base('eurusd','6MO','fx_eurusd_6mo','EUR/USD 6M')
fx_base('usdhkd','6MO','fx_usdhkd_6mo','USD/HKD 6M')

# eta sensitivity: EUR/USD 6M
for mode,mat,name,title in [('eurusd','6MO','fx_eurusd_eta_compare','EUR/USD 6M')]:
    s10=fxdata(FX10,mode,'continuous_smiles');s02=fxdata(FX02,mode,'continuous_smiles');d10=fxdata(FX10,mode,'densities');d02=fxdata(FX02,mode,'densities');q=fxdata(FX10,mode,'market_vs_models')
    s10=s10[s10.maturity.astype(str)==mat].sort_values('k');s02=s02[s02.maturity.astype(str)==mat].sort_values('k');d10=d10[d10.maturity.astype(str)==mat].sort_values('k');d02=d02[d02.maturity.astype(str)==mat].sort_values('k');q=q[q.maturity.astype(str)==mat].sort_values('k')
    klo=float(q.k.min());khi=float(q.k.max());w=khi-klo;xmin=max(float(s10.k.min()),klo-.85*w);xmax=min(float(s10.k.max()),khi+1.25*w)
    fig,axs=plt.subplots(1,2,figsize=(10.6,4.0))
    ax=axs[0];ax.plot(s10.k,100*s10.model1_lnmix_iv,lw=1.9,label=r'Mixture $\eta=10\%$');ax.plot(s02.k,100*s02.model1_lnmix_iv,lw=1.8,ls='-.',label=r'Mixture $\eta=2\%$');ax.plot(s10.k,100*s10.model2_cubic_iv,lw=1.4,ls='--',label='Cubic spline');ax.scatter(q.k,100*q.market_iv,s=11,marker='s',zorder=5,label='Market quotes');ax.axvline(klo,ls=':',lw=.8);ax.axvline(khi,ls=':',lw=.8);ax.set_xlim(xmin,xmax);ax.set_xlabel('Moneyness $K/F$');ax.set_ylabel('Implied volatility (%)');ax.set_title(title+' - floor sensitivity');ax.grid(alpha=.2);ax.legend()
    ax=axs[1];ax.plot(d10.k,d10.model1_density,lw=1.9,label=r'Density $\eta=10\%$');ax.plot(d02.k,d02.model1_density,lw=1.8,ls='-.',label=r'Density $\eta=2\%$');ax.plot(d10.k,d10.model2_density,lw=1.4,ls='--',label='Cubic-spline density');ax.axhline(0,lw=.7);ax.axvline(klo,ls=':',lw=.8);ax.axvline(khi,ls=':',lw=.8);ax.set_xlim(xmin,xmax);ax.set_xlabel('Moneyness $K/F$');ax.set_ylabel('Normalized risk-neutral density');ax.set_title(title+' - density sensitivity');ax.grid(alpha=.2);ax.legend();save(fig,name)

# USD/HKD all maturities: smiles and densities for eta 10% vs 2%
mats=['1MO','3MO','6MO','1YR']; labels={'1MO':'1M','3MO':'3M','6MO':'6M','1YR':'1Y'}
s10=fxdata(FX10,'usdhkd','continuous_smiles');s02=fxdata(FX02,'usdhkd','continuous_smiles');d10=fxdata(FX10,'usdhkd','densities');d02=fxdata(FX02,'usdhkd','densities');ma=fxdata(FX10,'usdhkd','market_vs_models')
fig,axs=plt.subplots(2,2,figsize=(10.5,7.2))
for ax,mat in zip(axs.ravel(),mats):
    a=s10[s10.maturity.astype(str)==mat].sort_values('k');b=s02[s02.maturity.astype(str)==mat].sort_values('k');q=ma[ma.maturity.astype(str)==mat].sort_values('k');klo=float(q.k.min());khi=float(q.k.max());w=khi-klo;xmin=max(float(a.k.min()),klo-.80*w);xmax=min(float(a.k.max()),khi+1.20*w)
    ax.plot(a.k,100*a.model1_lnmix_iv,lw=1.8,label=r'$\eta=10\%$');ax.plot(b.k,100*b.model1_lnmix_iv,lw=1.7,ls='-.',label=r'$\eta=2\%$');ax.plot(a.k,100*a.model2_cubic_iv,lw=1.25,ls='--',label='Cubic spline');ax.scatter(q.k,100*q.market_iv,s=10,marker='s',zorder=5,label='Market');ax.axvline(klo,ls=':',lw=.7);ax.axvline(khi,ls=':',lw=.7);ax.set_xlim(xmin,xmax);ax.set_title('USD/HKD '+labels[mat]);ax.set_xlabel('Moneyness $K/F$');ax.set_ylabel('IV (%)');ax.grid(alpha=.18)
handles,legs=axs[0,0].get_legend_handles_labels();fig.legend(handles,legs,loc='lower center',ncol=4,bbox_to_anchor=(.5,-.005));fig.suptitle('USD/HKD implied-volatility smiles - all quoted maturities',y=.995);fig.tight_layout(rect=[0,.045,1,.97]);fig.savefig(OUT/'fx_usdhkd_all_maturities_eta.pdf',bbox_inches='tight');fig.savefig(OUT/'fx_usdhkd_all_maturities_eta.png',dpi=190,bbox_inches='tight');plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(10.5,7.2))
for ax,mat in zip(axs.ravel(),mats):
    a=d10[d10.maturity.astype(str)==mat].sort_values('k');b=d02[d02.maturity.astype(str)==mat].sort_values('k');q=ma[ma.maturity.astype(str)==mat].sort_values('k');klo=float(q.k.min());khi=float(q.k.max());w=khi-klo;xmin=max(float(a.k.min()),klo-.80*w);xmax=min(float(a.k.max()),khi+1.20*w)
    ax.plot(a.k,a.model1_density,lw=1.8,label=r'Mixture $\eta=10\%$');ax.plot(b.k,b.model1_density,lw=1.7,ls='-.',label=r'Mixture $\eta=2\%$');ax.plot(a.k,a.model2_density,lw=1.25,ls='--',label='Cubic spline');ax.axhline(0,lw=.7);ax.axvline(klo,ls=':',lw=.7);ax.axvline(khi,ls=':',lw=.7);ax.set_xlim(xmin,xmax);ax.set_title('USD/HKD '+labels[mat]);ax.set_xlabel('Moneyness $K/F$');ax.set_ylabel('Density');ax.grid(alpha=.18)
handles,legs=axs[0,0].get_legend_handles_labels();fig.legend(handles,legs,loc='lower center',ncol=3,bbox_to_anchor=(.5,-.005));fig.suptitle('USD/HKD risk-neutral densities - all quoted maturities',y=.995);fig.tight_layout(rect=[0,.045,1,.97]);fig.savefig(OUT/'fx_usdhkd_all_maturities_density_eta.pdf',bbox_inches='tight');fig.savefig(OUT/'fx_usdhkd_all_maturities_density_eta.png',dpi=190,bbox_inches='tight');plt.close(fig)

# Equity baseline figures and eta comparison
def eqfiles(root,mode):
    return (pd.read_csv(root/mode/f'{mode}_continuous_smiles.csv'),pd.read_csv(root/mode/f'{mode}_densities.csv'),pd.read_csv(root/mode/f'{mode}_market_vs_models.csv'))

def eq_base(mode,mat,name,title):
    s,d,q=eqfiles(EQ10,mode);s=s[s.maturity.astype(str)==mat].sort_values('k');d=d[d.maturity.astype(str)==mat].sort_values('k');q=q[q.maturity.astype(str)==mat].sort_values('k');klo=float(q.k.min());khi=float(q.k.max());w=khi-klo;xmin=max(float(s.k.min()),klo-.08*w);xmax=min(float(s.k.max()),khi+.12*w)
    fig,axs=plt.subplots(1,2,figsize=(10.6,4.0));ax=axs[0];ax.plot(s.k,100*s.model1_flat_iv,lw=2,label=r'Flat mixture ($\eta=10\%$)');ax.plot(s.k,100*s.model2_cubic_iv,lw=1.6,ls='--',label='Cubic total-variance spline');ax.scatter(q.k,100*q.market_iv,s=9,marker='s',zorder=5,label='Market / benchmark knots');ax.axvline(klo,ls=':',lw=.8);ax.axvline(khi,ls=':',lw=.8);ax.set_xlim(xmin,xmax);ax.set_xlabel('Moneyness $K/F$');ax.set_ylabel('IV (%)');ax.set_title(title+' - market-region fit');ax.grid(alpha=.2);ax.legend()
    ax=axs[1];ax.plot(d.k,d.model1_flat_density,lw=2,label=r'Mixture density ($\eta=10\%$)');ax.plot(d.k,d.model2_density,lw=1.6,ls='--',label='Cubic-spline density');ax.axhline(0,lw=.7);ax.axvline(klo,ls=':',lw=.8);ax.axvline(khi,ls=':',lw=.8);ax.set_xlim(xmin,xmax);ax.set_xlabel('Moneyness $K/F$');ax.set_ylabel('Density');ax.set_title(title+' - density');ax.grid(alpha=.2);ax.legend();save(fig,name)

eq_base('spy_rich','2026-10-16','equity_spy_flat','SPY')
eq_base('tsla_rich','2026-10-16','equity_tsla_flat','TSLA')
eq_base('eurostoxx50','2014-10-17','equity_eurostoxx_flat','EURO STOXX 50')

cases=[('spy_rich','2026-10-16','SPY'),('tsla_rich','2026-10-16','TSLA'),('nvda_rich','2026-10-16','NVDA'),('aapl_rich','2026-10-16','AAPL'),('eurostoxx50','2014-10-17','EURO STOXX 50')]
fig,axs=plt.subplots(3,2,figsize=(10.6,10.3));axs=axs.ravel()
for ax,(mode,mat,title) in zip(axs,cases):
    s10,_,q=eqfiles(EQ10,mode);s02,_,_=eqfiles(EQ02,mode);s10=s10[s10.maturity.astype(str)==mat].sort_values('k');s02=s02[s02.maturity.astype(str)==mat].sort_values('k');q=q[q.maturity.astype(str)==mat].sort_values('k');klo=float(q.k.min());khi=float(q.k.max());w=khi-klo;xmin=klo-.04*w;xmax=khi+.06*w
    ax.plot(s10.k,100*s10.model1_flat_iv,lw=1.7,label=r'$\eta=10\%$');ax.plot(s02.k,100*s02.model1_flat_iv,lw=1.7,ls='-.',label=r'$\eta=2\%$');ax.scatter(q.k,100*q.market_iv,s=8,marker='s',zorder=5,label='Market');ax.set_xlim(xmin,xmax);ax.set_title(title);ax.set_xlabel('$K/F$');ax.set_ylabel('IV (%)');ax.grid(alpha=.18)
axs[-1].axis('off');handles,legs=axs[0].get_legend_handles_labels();fig.legend(handles,legs,loc='lower center',ncol=3,bbox_to_anchor=(.5,.01));fig.suptitle('Equity/index market-region fit: sensitivity to the kernel floor',y=.995);fig.tight_layout(rect=[0,.04,1,.97]);fig.savefig(OUT/'equity_eta_smiles.pdf',bbox_inches='tight');fig.savefig(OUT/'equity_eta_smiles.png',dpi=190,bbox_inches='tight');plt.close(fig)

fig,axs=plt.subplots(3,2,figsize=(10.6,10.3));axs=axs.ravel()
for ax,(mode,mat,title) in zip(axs,cases):
    _,d10,q=eqfiles(EQ10,mode);_,d02,_=eqfiles(EQ02,mode);d10=d10[d10.maturity.astype(str)==mat].sort_values('k');d02=d02[d02.maturity.astype(str)==mat].sort_values('k');q=q[q.maturity.astype(str)==mat].sort_values('k');klo=float(q.k.min());khi=float(q.k.max());w=khi-klo;xmin=klo-.04*w;xmax=khi+.06*w
    ax.plot(d10.k,d10.model1_flat_density,lw=1.7,label=r'$\eta=10\%$');ax.plot(d02.k,d02.model1_flat_density,lw=1.7,ls='-.',label=r'$\eta=2\%$');ax.axhline(0,lw=.7);ax.set_xlim(xmin,xmax);ax.set_title(title);ax.set_xlabel('$K/F$');ax.set_ylabel('Density');ax.grid(alpha=.18)
axs[-1].axis('off');handles,legs=axs[0].get_legend_handles_labels();fig.legend(handles,legs,loc='lower center',ncol=2,bbox_to_anchor=(.5,.01));fig.suptitle('Equity/index mixture density: sensitivity to the kernel floor',y=.995);fig.tight_layout(rect=[0,.04,1,.97]);fig.savefig(OUT/'equity_eta_densities.pdf',bbox_inches='tight');fig.savefig(OUT/'equity_eta_densities.png',dpi=190,bbox_inches='tight');plt.close(fig)

# stress plots, with small square markers and multiple SPX alpha values
stress_market=pd.read_csv(STRESS/'all_market_vs_models.csv');stress_wings=pd.read_csv(STRESS/'all_wing_extrapolation.csv');comps=pd.read_csv(STRESS/'all_model1_components.csv')
def ncdf(x): return .5*(1+math.erf(x/math.sqrt(2)))
def sc(center,k,v):
    if v<=1e-12:return max(center-k,0)
    d1=(-math.log(k/center)+.5*v*v)/v;d2=d1-v;return center*ncdf(d1)-k*ncdf(d2)
def pc(k,v):
    if v<=1e-12:return max(1-k,0)
    d1=(-math.log(k)+.5*v*v)/v;d2=d1-v;return ncdf(d1)-k*ncdf(d2)
def imp(k,c,t):
    lo,hi=1e-9,3
    target=max(c,max(1-k,0)+1e-14)
    for _ in range(100):
        mid=.5*(lo+hi)
        if pc(k,mid)<target:lo=mid
        else:hi=mid
    return .5*(lo+hi)/math.sqrt(t)
def alpha_smile(inst,mat,kgrid,frac):
    c=comps[(comps.instrument==inst)&(comps.maturity.astype(str)==mat)].copy();alpha=float(c.right_wing_alpha.iloc[0])*frac;atm=float(c.atm_total_vol.iloc[0]);kmax=float(c.market_k_max.iloc[0]);t=float(c.ttm.iloc[0]);out=[]
    for k in kgrid:
        call=0
        for _,r in c.iterrows():
            base=float(r.kernel_total_vol);sig=base
            if k>kmax:sig=base*math.exp(-alpha*(base/max(atm,1e-14))**2*math.log(k/kmax)**2)
            call+=float(r.weight)*sc(float(r.center_k),float(k),sig)
        out.append(imp(float(k),call,t))
    return np.array(out)
def stress(inst,name,title):
    q=stress_market[stress_market.instrument==inst].sort_values('k');mat=str(q.maturity.iloc[0]);g=stress_wings[(stress_wings.instrument==inst)&(stress_wings.maturity.astype(str)==mat)].sort_values('k');klo=float(q.k.min());khi=float(q.k.max());w=khi-klo;xmin=max(float(g.k.min()),klo-.10*w);xmax=float(g.k.max())
    fig,axs=plt.subplots(1,2,figsize=(10.6,4.0));ax=axs[0];ax.plot(g.k,100*g.model1_flat_iv,lw=1.7,ls='--',label='Flat widths')
    if inst=='SPX':
        kg=g.k.to_numpy(float);ax.plot(kg,100*alpha_smile(inst,mat,kg,.25),lw=1.7,label=r'$0.25\alpha_R$');ax.plot(kg,100*alpha_smile(inst,mat,kg,.50),lw=1.7,label=r'$0.50\alpha_R$');ax.plot(g.k,100*g.model1_winged_iv,lw=1.8,label=r'Full $\alpha_R$')
    else:ax.plot(g.k,100*g.model1_winged_iv,lw=1.9,label='Controlled right wing')
    ax.scatter(q.k,100*q.market_iv,s=9,marker='s',zorder=5,label='Market');ax.axvline(khi,ls=':',lw=.8,label='Last market strike');ax.set_xlim(xmin,xmax);ax.set_xlabel('$K/F$');ax.set_ylabel('IV (%)');ax.set_title(title+' - right wing');ax.grid(alpha=.2);ax.legend()
    ax=axs[1];ax.plot(g.k,g.model1_flat_density,lw=1.7,ls='--',label='Flat density');ax.plot(g.k,g.model1_winged_density,lw=1.9,label='Full-control density');ax.axhline(0,lw=.7);ax.axvline(khi,ls=':',lw=.8);ax.set_xlim(xmin,xmax);ax.set_xlabel('$K/F$');ax.set_ylabel('Density');ax.set_title(title+' - extrapolated density');ax.grid(alpha=.2);ax.legend();save(fig,name)
stress('SPX','wing_spx_stress','SPX');stress('TSLA','wing_tsla_stress','TSLA')
print('Figures written to',OUT)
