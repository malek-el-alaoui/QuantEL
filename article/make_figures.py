from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

OUT = Path('/mnt/data/revised_lognormal_mixture_article_v3/figures')
FXROOT = Path('/mnt/data/fx_flat_with_usdhkd/output')
EQROOT = Path('/mnt/data/VolatilitySmileTwoModels_RichEurostox30_BetterPlots/output')

plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'legend.fontsize': 8.5,
    'figure.dpi': 140,
})

def bounds_from_market(klo,khi,left=0.20,right=0.28):
    w=max(khi-klo,1e-6)
    return max(1e-6,klo-left*w), khi+right*w

def save_fx(prefix, maturity, name, title):
    sm=pd.read_csv(FXROOT/f'{prefix}_continuous_smiles.csv')
    de=pd.read_csv(FXROOT/f'{prefix}_densities.csv')
    ma=pd.read_csv(FXROOT/f'{prefix}_market_vs_models.csv')
    s=sm[sm.maturity.astype(str)==str(maturity)].sort_values('k')
    d=de[de.maturity.astype(str)==str(maturity)].sort_values('k')
    q=ma[ma.maturity.astype(str)==str(maturity)].sort_values('k')
    klo=float(q.k.min()); khi=float(q.k.max()); w=khi-klo
    # use essentially full 0.5-delta model grid, but slightly more room on right
    xmin=max(float(s.k.min()), klo-1.05*w)
    xmax=min(float(s.k.max()), khi+1.45*w)
    fig,axs=plt.subplots(1,2,figsize=(10.6,4.15))
    ax=axs[0]
    ax.plot(s.k,100*s.model1_lnmix_iv,lw=2,label='Log-normal mixture')
    ax.plot(s.k,100*s.model2_cubic_iv,lw=1.7,ls='--',label='Cubic total-variance spline')
    ax.scatter(q.k,100*q.market_iv,s=20,marker='s',zorder=5,label='Market quotes')
    ax.axvline(klo,ls=':',lw=0.9); ax.axvline(khi,ls=':',lw=0.9)
    ax.set_xlim(xmin,xmax); ax.set_xlabel('Moneyness $K/F$'); ax.set_ylabel('Implied volatility (%)')
    ax.set_title(title+' - smile'); ax.grid(alpha=.22); ax.legend()
    ax=axs[1]
    ax.plot(d.k,d.model1_density,lw=2,label='Log-normal-mixture density')
    ax.plot(d.k,d.model2_density,lw=1.7,ls='--',label='Cubic-spline density')
    ax.axhline(0,lw=.8); ax.axvline(klo,ls=':',lw=.9); ax.axvline(khi,ls=':',lw=.9)
    ax.set_xlim(xmin,xmax); ax.set_xlabel('Moneyness $K/F$'); ax.set_ylabel('Normalized risk-neutral density')
    ax.set_title(title+' - density'); ax.grid(alpha=.22); ax.legend()
    fig.tight_layout()
    fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight')
    fig.savefig(OUT/f'{name}.png',dpi=180,bbox_inches='tight')
    plt.close(fig)

# FX figures
save_fx('eurusd','6MO','fx_eurusd_6mo','EUR/USD 6M')
save_fx('usdhkd','6MO','fx_usdhkd_6mo','USD/HKD 6M')

# equity base figures: FLAT mixture only in prior results section
market=pd.read_csv(EQROOT/'all_market_vs_models.csv')
smiles=pd.read_csv(EQROOT/'all_continuous_smiles.csv')
dens=pd.read_csv(EQROOT/'all_densities.csv')
wings=pd.read_csv(EQROOT/'all_wing_extrapolation.csv')

def equity_base(inst,mat,name,title):
    q=market[(market.instrument==inst)&(market.maturity.astype(str)==str(mat))].sort_values('k')
    s=smiles[(smiles.instrument==inst)&(smiles.maturity.astype(str)==str(mat))].sort_values('k')
    d=dens[(dens.instrument==inst)&(dens.maturity.astype(str)==str(mat))].sort_values('k')
    klo=float(q.k.min()); khi=float(q.k.max()); w=khi-klo
    xmin=max(float(s.k.min()),klo-.10*w); xmax=min(float(s.k.max()),khi+.16*w)
    fig,axs=plt.subplots(1,2,figsize=(10.6,4.15))
    ax=axs[0]
    ax.plot(s.k,100*s.model1_flat_iv,lw=2,label='Flat-kernel mixture')
    ax.plot(s.k,100*s.model2_cubic_iv,lw=1.7,ls='--',label='Cubic total-variance spline')
    ax.scatter(q.k,100*q.market_iv,s=18,marker='s',zorder=5,label='Market / benchmark knots')
    ax.axvline(klo,ls=':',lw=.9); ax.axvline(khi,ls=':',lw=.9)
    ax.set_xlim(xmin,xmax); ax.set_xlabel('Moneyness $K/F$'); ax.set_ylabel('Implied volatility (%)')
    ax.set_title(title+' - market-region fit'); ax.grid(alpha=.22); ax.legend()
    ax=axs[1]
    ax.plot(d.k,d.model1_flat_density,lw=2,label='Flat-kernel mixture density')
    ax.plot(d.k,d.model2_density,lw=1.7,ls='--',label='Cubic-spline density')
    ax.axhline(0,lw=.8); ax.axvline(klo,ls=':',lw=.9); ax.axvline(khi,ls=':',lw=.9)
    ax.set_xlim(xmin,xmax); ax.set_xlabel('Moneyness $K/F$'); ax.set_ylabel('Normalized risk-neutral density')
    ax.set_title(title+' - density'); ax.grid(alpha=.22); ax.legend()
    fig.tight_layout()
    fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight')
    fig.savefig(OUT/f'{name}.png',dpi=180,bbox_inches='tight')
    plt.close(fig)

equity_base('SPY_RICH','2026-10-16','equity_spy_flat','SPY')
equity_base('TSLA_RICH','2026-10-16','equity_tsla_flat','TSLA')
equity_base('EUROSTOXX50','2014-10-17','equity_eurostoxx_flat','EURO STOXX 50')

# separate wing-extrapolation figures, show flat vs controlled only + market points

def equity_wing(inst,mat,name,title):
    q=market[(market.instrument==inst)&(market.maturity.astype(str)==str(mat))].sort_values('k')
    g=wings[(wings.instrument==inst)&(wings.maturity.astype(str)==str(mat))].sort_values('k')
    klo=float(q.k.min()); khi=float(q.k.max()); w=khi-klo
    # deliberately suppress most of left extrapolation and display substantially more right wing
    xmin=max(float(g.k.min()), klo-.08*w)
    xmax=min(float(g.k.max()), khi+.36*w)
    fig,axs=plt.subplots(1,2,figsize=(10.6,4.15))
    ax=axs[0]
    ax.plot(g.k,100*g.model1_flat_iv,lw=1.8,ls='--',label='Flat kernel widths')
    ax.plot(g.k,100*g.model1_winged_iv,lw=2.2,label='Controlled equity right wing')
    ax.scatter(q.k,100*q.market_iv,s=18,marker='s',zorder=5,label='Market quotes')
    ax.axvline(khi,ls=':',lw=1,label='Last market strike')
    ax.set_xlim(xmin,xmax); ax.set_xlabel('Moneyness $K/F$'); ax.set_ylabel('Implied volatility (%)')
    ax.set_title(title+' - right-wing extrapolation'); ax.grid(alpha=.22); ax.legend()
    ax=axs[1]
    ax.plot(g.k,g.model1_flat_density,lw=1.8,ls='--',label='Flat-kernel density')
    ax.plot(g.k,g.model1_winged_density,lw=2.2,label='Controlled-wing density')
    ax.axhline(0,lw=.8); ax.axvline(khi,ls=':',lw=1,label='Last market strike')
    ax.set_xlim(xmin,xmax); ax.set_xlabel('Moneyness $K/F$'); ax.set_ylabel('Normalized risk-neutral density')
    ax.set_title(title+' - extrapolated density'); ax.grid(alpha=.22); ax.legend()
    fig.tight_layout()
    fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight')
    fig.savefig(OUT/f'{name}.png',dpi=180,bbox_inches='tight')
    plt.close(fig)

equity_wing('SPY_RICH','2026-10-16','wing_spy','SPY')
equity_wing('TSLA_RICH','2026-10-16','wing_tsla','TSLA')
equity_wing('AAPL_RICH','2026-10-16','wing_aapl','AAPL')
equity_wing('NVDA_RICH','2026-10-16','wing_nvda','NVDA')

# Table-ready numeric summaries
# FX summary
fxma=[]
for prefix in ['eurusd','usdhkd']:
    m=pd.read_csv(FXROOT/f'{prefix}_market_vs_models.csv')
    den=pd.read_csv(FXROOT/f'{prefix}_densities.csv')
    for mat in sorted(m.maturity.astype(str).unique()):
        mm=m[m.maturity.astype(str)==mat]
        dd=den[den.maturity.astype(str)==mat]
        fxma.append((prefix.upper(),mat,len(mm),float(mm.model1_iv_error_bp.abs().max()),float(dd.model1_density.min()),float(dd.model2_density.min())))
print('FXSUM')
for r in fxma: print(r)
print('WINGSUM')
for inst in ['SPY_RICH','TSLA_RICH','NVDA_RICH','AAPL_RICH']:
    g=wings[wings.instrument==inst].sort_values('k')
    q=market[market.instrument==inst]
    khi=float(q.k.max())
    right=g[g.k>=khi]
    last=right.iloc[-1]
    b=right.iloc[0]
    print(inst, 'flat_end',float(last.model1_flat_iv),'ctrl_end',float(last.model1_winged_iv),'alpha',float(last.right_wing_alpha),'minctrl',float(right.model1_winged_density.min()))
# Additional stress-case wing figures from the bounded-kernel equity project,
# chosen because the uncontrolled flat-kernel right wing is visibly unrealistic.
STRESS = Path('/mnt/data/VolatilitySmileTwoModels_EquityBoundedKernels_RightWing/output')
stress_market = pd.read_csv(STRESS/'all_market_vs_models.csv')
stress_wings = pd.read_csv(STRESS/'all_wing_extrapolation.csv')

def _norm_cdf(x):
    import math
    return 0.5*(1.0+math.erf(x/math.sqrt(2.0)))

def _shifted_call(center,k,total_vol):
    import math
    if total_vol <= 1e-12:
        return max(center-k,0.0)
    d1=(-math.log(k/center)+0.5*total_vol*total_vol)/total_vol
    d2=d1-total_vol
    return center*_norm_cdf(d1)-k*_norm_cdf(d2)

def _plain_call(k,total_vol):
    import math
    if total_vol <= 1e-12:
        return max(1.0-k,0.0)
    d1=(-math.log(k)+0.5*total_vol*total_vol)/total_vol
    d2=d1-total_vol
    return _norm_cdf(d1)-k*_norm_cdf(d2)

def _implied_ann_vol(k,call,t):
    import math
    target=max(call,max(1.0-k,0.0)+1e-14)
    lo,hi=1e-9,3.0
    for _ in range(120):
        mid=0.5*(lo+hi)
        if _plain_call(k,mid) < target: lo=mid
        else: hi=mid
    return 0.5*(lo+hi)/math.sqrt(t)

def _scaled_alpha_smile(inst,mat,kgrid,alpha_fraction):
    import math
    comps=pd.read_csv(STRESS/'all_model1_components.csv')
    c=comps[(comps.instrument==inst)&(comps.maturity.astype(str)==str(mat))].copy()
    alpha=float(c.right_wing_alpha.iloc[0])*alpha_fraction
    atm=float(c.atm_total_vol.iloc[0]); kmax=float(c.market_k_max.iloc[0]); t=float(c.ttm.iloc[0])
    vals=[]
    for k in kgrid:
        cp=0.0
        for _,r in c.iterrows():
            base=float(r.kernel_total_vol)
            sig=base
            if k>kmax:
                fac=(base/max(atm,1e-14))**2
                sig=base*math.exp(-alpha*fac*math.log(k/kmax)**2)
            cp += float(r.weight)*_shifted_call(float(r.center_k),float(k),sig)
        vals.append(_implied_ann_vol(float(k),cp,t))
    return np.asarray(vals)

def stress_wing(inst,name,title):
    q=stress_market[stress_market.instrument==inst].sort_values('k')
    mat=str(q.maturity.iloc[0])
    g=stress_wings[(stress_wings.instrument==inst)&(stress_wings.maturity.astype(str)==mat)].sort_values('k')
    klo=float(q.k.min()); khi=float(q.k.max()); w=khi-klo
    xmin=max(float(g.k.min()),klo-.18*w)
    xmax=float(g.k.max())
    fig,axs=plt.subplots(1,2,figsize=(10.6,4.15))
    ax=axs[0]
    ax.plot(g.k,100*g.model1_flat_iv,lw=1.8,ls='--',label='Flat kernel widths')
    if inst=='SPX':
        kgrid=g.k.to_numpy(float)
        iv25=_scaled_alpha_smile(inst,mat,kgrid,0.25)
        iv50=_scaled_alpha_smile(inst,mat,kgrid,0.50)
        ax.plot(kgrid,100*iv25,lw=1.8,label=r'Mild control: $0.25\,\alpha_R$')
        ax.plot(kgrid,100*iv50,lw=1.9,label=r'Medium control: $0.50\,\alpha_R$')
        ax.plot(g.k,100*g.model1_winged_iv,lw=2.0,label=r'Full control: $\alpha_R$')
    else:
        ax.plot(g.k,100*g.model1_winged_iv,lw=2.2,label='Controlled right wing')
    ax.scatter(q.k,100*q.market_iv,s=18,marker='s',zorder=5,label='Market quotes')
    ax.axvline(khi,ls=':',lw=1,label='Last market strike')
    ax.set_xlim(xmin,xmax); ax.set_xlabel('Moneyness $K/F$'); ax.set_ylabel('Implied volatility (%)')
    ax.set_title(title+' - right-wing extrapolation'); ax.grid(alpha=.22); ax.legend()
    ax=axs[1]
    ax.plot(g.k,g.model1_flat_density,lw=1.8,ls='--',label='Flat-kernel density')
    ax.plot(g.k,g.model1_winged_density,lw=2.2,label='Full-control density')
    ax.axhline(0,lw=.8); ax.axvline(khi,ls=':',lw=1,label='Last market strike')
    ax.set_xlim(xmin,xmax); ax.set_xlabel('Moneyness $K/F$'); ax.set_ylabel('Normalized risk-neutral density')
    ax.set_title(title+' - extrapolated density'); ax.grid(alpha=.22); ax.legend()
    fig.tight_layout()
    fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight')
    fig.savefig(OUT/f'{name}.png',dpi=180,bbox_inches='tight')
    plt.close(fig)
    last=g.iloc[-1]
    print('STRESS', inst, mat, 'flat_end',float(last.model1_flat_iv),'controlled_end',float(last.model1_winged_iv), 'alpha', float(last.right_wing_alpha))

stress_wing('SPX','wing_spx_stress','SPX')
stress_wing('TSLA','wing_tsla_stress','TSLA')
