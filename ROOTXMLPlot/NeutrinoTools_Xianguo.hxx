#include <fstream>

#include <TLorentzVector.h>
#include <TH3D.h>
#include <TProfile2D.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>

const Double_t EPSILON = 1e-12;

class NeutrinoTools
{
 public:
  enum errflag{
    kIniValue    = -999,
    kNullPointer = -9999,
    kZeroDivider = -99999
  };

  //==========================================================
  //General constants
  //==========================================================
  //http://www.star.bnl.gov/public/comp/simu/gstar/Manual/particle_id.html
  //static UInt_t GeantNeutrino(){ return 4;}
  //static UInt_t GeantMuonPlus(){ return 5;}
  //static UInt_t GeantMuonMinus(){ return 6;}
  //static UInt_t GeantPionPlus(){return 8;}
  //static UInt_t GeantPionMinus(){return 9;}
  //http://pdg.lbl.gov/2011/reviews/rpp2011-rev-monte-carlo-numbering.pdf
  static Int_t PDGNumu(){ return 14;}
  static Int_t PDGAntiNumu(){ return -14;}
  static Int_t PDGMuonPlus(){ return -13;}
  static Int_t PDGMuonMinus(){ return 13;}
  static Int_t PDGPionPlus(){ return 211;}
  static Int_t PDGPionMinus(){ return -211;}
  static Int_t PDGKLong(){ return 130;}
  static Int_t PDGKaonPlus(){ return 321;}
  static Int_t PDGKaonMinus(){ return -321;}

  static Double_t MuonMass(){ return 105.65837/1e3; }//in GeV //google = wiki
  static Double_t ProtonMass(){ return 938.272/1e3;}//in GeV //google = wiki
  static Double_t PionMass(){ return 139.570/1e3;}//in GeV //wiki
  static Double_t NeutronMass(){ return 939.565/1e3;}//in GeV //wiki
  static Double_t DeltaPPMass(){ return 1232/1e3;}//in GeV// resonance, only approximate

  //http://journals.aps.org/prd/pdf/10.1103/PhysRevD.88.032002
  //"and the binding energy is set to 25 MeV for carbon
  //and 27 MeV for oxygen"
  static Double_t BindingEnergyCarbon(){return 25/1e3;}//in GeV

  static Int_t PDGToType(const Int_t pdg);

  //==========================================================
  //General algorithms
  //==========================================================

  static Double_t GetRoughMPVWithFWHM(TH1 * hh, Double_t & lowerX, Double_t & higherX);
  static TH1D*    ToPDF(const TH1 *hraw, const TString hn="");
  static TH2*    NormalHist(const TH2 *hraw, const Double_t thres=0, const Bool_t kmax=kFALSE);
  static void     ToNaturalScale(TAxis *ax);
  static void     ToNaturalScale(TH1 *hh);
  static void     BinLog(TAxis *axis, const Double_t non0start=-999);
  static TH2D *   ProjectionYX(const TH3D *hh, const Bool_t klogx, const Bool_t klogy, const Int_t iz0, const Int_t iz1, Int_t &count);
  static TGraphAsymmErrors * GetFluxMap(const Double_t tol, const TH2D *hh, TH2D * & hmap, TH2D * & hdiff);
  static TChain * InputFiles(const TString file, const TString tr, Char_t *dir=0x0);
  static Double_t * GetAxisArray(TAxis * aa);
  static void FitSlicesY(const TH2D *hh, TH1D *&hnor, TH1D *&hmpv, TH1D *&hwid, TH1D *&hres, TH1D *&hchi, const TString formula, const Double_t thres, TList *ll=0x0);
  static void ScaleToRef(TH1D * hh, const TH1D *href);

  //==========================================================
  //General calculations
  //==========================================================
  static Double_t SampleFermiMomentum();
  static void SampleUnityIsotropicVector(TVector3 * vec);
  static void Deflection(TVector3 * vec, const Double_t dtheta, const Double_t dphi);
  static Double_t RutherfordTheta(const Double_t regpar);
  static Double_t RutherfordEnergy(const Double_t regpar);
  static void RutherfordTransport(TLorentzVector * v0, const Int_t nstep, const Double_t regphi, const Double_t regtheta, const Double_t regde);

  static Double_t GetPtFast(const TVector3 * refdir, const Double_t dir0, const Double_t dir1, const Double_t dir2, const Double_t pp);
  static const TVector3 *GetVecT(const TLorentzVector * refdir, const Double_t xx, const Double_t yy, const Double_t zz);
  static Double_t CalcAlpha(const TLorentzVector *refdir, const TLorentzVector * mup4);
  static Double_t GetCos(const TVector3 * v1, const TVector3 * v2, const TString tag);
  static Double_t GetCos(const TLorentzVector * v1, const TLorentzVector * v2, const TString tag){ const TVector3 tmp1=v1->Vect(); const TVector3 tmp2=v2->Vect();  return GetCos(&tmp1, &tmp2, tag+"LorentzVector"); }
  static Double_t GetAngle(const TVector3 * v1, const TVector3 * v2, const TString tag);
  static Double_t GetAngle(const TLorentzVector * v1, const TLorentzVector * v2, const TString tag){ const TVector3 tmp1=v1->Vect(); const TVector3 tmp2=v2->Vect();  return GetAngle(&tmp1, &tmp2, tag+"LorentzVector"); }
  static Double_t GetSin(const TVector3 * v1, const TVector3 * v2, const TString tag);
  static Double_t GetSin(const TLorentzVector * v1, const TLorentzVector * v2, const TString tag){ const TVector3 tmp1=v1->Vect(); const TVector3 tmp2=v2->Vect();  return GetSin(&tmp1, &tmp2, tag+"LorentzVector"); }
  static Double_t GetDeltaPTT(const TVector3 &p1, const TVector3 &p2, const TVector3 & nuDir, const TVector3 & refDir);
  static Double_t GetMass(const TLorentzVector * par1, const TLorentzVector * par2);
  static Double_t GetPTT(const TVector3 & parMom, const TVector3 & nuDir, const TVector3 & muonDir, const Bool_t kprint=kFALSE);
  //static Double_t GetKT(const TVector3 * ptScal, const TVector3 * ptAxis);
  static void SetDeltaPt(TVector3 * deltapt, const TVector3 * ptmu, const TVector3 * ptproton);

  static Double_t GetEnuApp(const TVector3 & v3Nu, const Double_t mN, const TLorentzVector *lvL, const Double_t mX);

  //==========================================================
  //IO
  //==========================================================
  static Int_t FillTH1I(TList *lin, const Int_t id, const Int_t var);
  static Int_t FillTH1D(TList *lin, const Int_t id, const Double_t var);
  static Int_t FillTH2D(TList *lin, const Int_t id, const Double_t var1, const Double_t var2);
  static Int_t FillTH3D(TList *lin, const Int_t id, const Double_t var1, const Double_t var2, const Double_t var3);

  //==========================================================
  //Plotting related
  //==========================================================
  static TString GetTitleFromVar(const TString var);
  //static TString GetPIDCut();
  //static TString GetPtResCut();
  //static TString GetQACut();
  static TString GetModeCut(const Int_t imode, TString & name);

 private:
  static Double_t GetCosDeltaPhiT(const TVector3 * ptmuon, const TVector3 * ptproton);

  static TRandom3 * fRan;
};

TRandom3 *NeutrinoTools::fRan=new TRandom3(123);

Int_t NeutrinoTools::PDGToType(const Int_t pdg)
{
  Int_t type = kIniValue;

  if(abs(pdg)==11){
    type = 1;
  }
  else if(abs(pdg)==13){
    type = 2;
  }
  else if(abs(pdg)==130){
    type = 3;
  }
  else if(abs(pdg)==211){
    type = 4;
  }
  else if(abs(pdg)==321){
    type = 5;
  }
  else if(abs(pdg)==2212){
    type = 6;
  }
  else if(abs(pdg)!=999 && pdg<1E8){
    printf("NeutrinoTools::PDGToType bad pdg %d\n", pdg); exit(1);
  }

  return type*(pdg>0?1:-1);
}

Double_t NeutrinoTools::GetRoughMPVWithFWHM(TH1 * hh, Double_t & lowerX, Double_t & higherX)
{
  const Double_t hmax = hh->GetBinContent(hh->GetMaximumBin());
  const Int_t bin1 = hh->FindFirstBinAbove(hmax/2.);
  const Int_t bin2 = hh->FindLastBinAbove(hmax/2.);
  lowerX = hh->GetBinCenter(bin1);
  higherX = hh->GetBinCenter(bin2);
  return hh->GetBinCenter(hh->GetMaximumBin());
}

TH1D* NeutrinoTools::ToPDF(const TH1 *hraw, const TString hn)
{
  const Int_t x0 = 0;
  const Int_t x1 = hraw->GetNbinsX()+1;
  const Double_t tmpnt = hraw->Integral(x0, x1);

  TH1D * hist = (TH1D*) hraw->Clone((hn+hraw->GetName())+"pdf");
  hist->Scale(0);

  /*
  //only consider region in range, ignore over/underflow. Important for fitting experiment data where under/overflow is ignored.
  const Int_t x0 = hist->GetXaxis()->GetFirst();
  const Int_t x1 = hist->GetXaxis()->GetLast();

  for(Int_t ii=x0; ii<=x1; ii++){
    const Double_t err = hraw->GetBinError(ii);
    const Double_t cont = hraw->GetBinContent(ii);

    //skip empty bins
    if(err<EPSILON){
      if(cont>EPSILON){
        printf("NeutrinoTools::ToPDF error! %d %e %e\n", ii, err, cont); exit(1);
      }
      continue;
    }

    const Double_t nn = cont*cont/err/err;
    hist->SetBinContent(ii,nn);
    hist->SetBinError(ii, TMath::Sqrt(nn));
    tmpnt += nn;
  }

  if(tmpnt<EPSILON){
    printf("NeutrinoTools::ToPDF tmpnt<epsilon ! %f %f\n", hist->GetEntries(), tmpnt);
    exit(1);
  }
  */

  for(Int_t ib=x0; ib<=x1; ib++){
    const Double_t bw = hraw->GetBinWidth(ib);
    const Double_t cont = hraw->GetBinContent(ib);
    if(cont<EPSILON)
      continue;

    //in case of finit number of bins (i.e. eff not always small), Binomial error is more accurate than Poisson error
    const Double_t eff = cont/tmpnt;
    const Double_t pdf = eff/bw;

    const Double_t dpdf = sqrt(eff*(1-eff)/tmpnt) / bw;
    hist->SetBinContent(ib, pdf);
    hist->SetBinError(ib, dpdf);
  }

  hist->SetEntries(tmpnt);

  return hist;
}

//thresh, min integral to compute
//kmax true: each slice the max is at 1,
TH2* NeutrinoTools::NormalHist(const TH2 *hraw, const Double_t thres, const Bool_t kmax) {

  TH2 *hh=(TH2*)hraw->Clone(Form("%snor",hraw->GetName()));
  hh->Scale(0);

  const Int_t x0 = hh->GetXaxis()->GetFirst();
  const Int_t x1 = hh->GetXaxis()->GetLast();
  const Int_t y0 = hh->GetYaxis()->GetFirst();
  const Int_t y1 = hh->GetYaxis()->GetLast();

  Double_t hmax = -1e10;
  Double_t hmin = 1e10;
  Double_t nent = 0;
  for(Int_t ix=x0; ix<=x1; ix++){

    //if option "e" is specified, the errors are computed. if option "o" original axis range of the taget axes will be kept, but only bins inside the selected range will be filled.

 TH1 * sliceh = hraw->ProjectionY(Form("tmpnormalhist%sx%d", hh->GetName(), ix), ix, ix, "oe");
 const Double_t tot = sliceh->GetEntries();

 TH1 * pdfh=0x0;

 if(tot>EPSILON){
   nent += tot;

   Double_t imax = kIniValue;

   if(!kmax){
     pdfh = ToPDF(sliceh,"tmp");
   }
   else{
     imax = sliceh->GetBinContent(sliceh->GetMaximumBin());
   }

   for(Int_t iy=y0; iy<=y1; iy++){
     const Double_t cont = kmax ? sliceh->GetBinContent(iy)/imax : pdfh->GetBinContent(iy);
     const Double_t ierr = kmax ? sliceh->GetBinError(iy)/imax   : pdfh->GetBinError(iy);
     if(tot>thres && cont>0){
       hh->SetBinContent(ix, iy, cont);
       hh->SetBinError(ix,iy, ierr);
       if(cont>hmax) hmax = cont;
       if(cont<hmin) hmin = cont;
     }
   }
 }

 delete pdfh;
 delete sliceh;
  }

  hh->SetEntries(nent);
  hh->SetMinimum(0.99*hmin);
  hh->SetMaximum(1.1*hmax);

  TString xtit(hraw->GetXaxis()->GetTitle());
  if(xtit.Contains("(")){
    xtit=xtit(0, xtit.First('('));
  }

  TString ytit(hraw->GetYaxis()->GetTitle());
  if(ytit.Contains("(")){
    ytit=ytit(0, ytit.First('('));
  }

  hh->SetTitle(Form("f(%s|%s) %s", ytit.Data(), xtit.Data(), hraw->GetTitle()));
  return hh;
}

void NeutrinoTools::ToNaturalScale(TAxis *ax)
{
  TAxis* oldx = (TAxis*)ax->Clone("oldx");
  ax->SetLimits(TMath::Power(10,oldx->GetXmin()), TMath::Power(10,oldx->GetXmax()));
  const Int_t nb = oldx->GetNbins();
  Double_t *bins = new Double_t[nb+1];
  bins[0]=TMath::Power(10,oldx->GetXmin());
  for(Int_t ii=1; ii<=nb; ii++){
    bins[ii]=TMath::Power(10,oldx->GetBinUpEdge(ii));
  }
  ax->Set(nb, bins);

  delete oldx;
  delete bins;
}

void NeutrinoTools::ToNaturalScale(TH1 *hh)
{
  ToNaturalScale(hh->GetXaxis());
}

void NeutrinoTools::BinLog(TAxis *axis, const Double_t non0start)
{
  const Int_t bins = axis->GetNbins();

  const Double_t xmin = axis->GetXmin();
  const Double_t xmax = axis->GetXmax();

  Bool_t k0start = kFALSE;
  if (xmin<EPSILON){
    k0start = kTRUE;
    if(non0start<EPSILON){
      printf("NeutrinoTools::BinLog bad non0start %f\n", non0start); exit(1);
    }
  }

  Double_t *new_bins = new Double_t[bins + 1];

  const Double_t factor = k0start? (TMath::Power(xmax/non0start, 1./(bins-1))) : (TMath::Power(xmax/xmin, 1./bins)) ;

  new_bins[0] = xmin;
  new_bins[1] = k0start ? non0start : (new_bins[0]*factor);

  for (int i = 2; i <= bins; i++) {
    new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete [] new_bins;
}

TH2D * NeutrinoTools::ProjectionYX(const TH3D *hh, const Bool_t klogx, const Bool_t klogy, const Int_t iz0, const Int_t iz1, Int_t &count)
{
  TH2D *h2= new TH2D(Form("%sproj%d", hh->GetName(), count++), "",
                     hh->GetNbinsX(), hh->GetXaxis()->GetBinLowEdge(1), hh->GetXaxis()->GetBinUpEdge(hh->GetNbinsX()),
                     hh->GetNbinsY(), hh->GetYaxis()->GetBinLowEdge(1), hh->GetYaxis()->GetBinUpEdge(hh->GetNbinsY()) );
  if(klogx){
    BinLog(h2->GetXaxis());
  }
  if(klogy){
    BinLog(h2->GetYaxis());
  }

  TString ztit(hh->GetZaxis()->GetTitle());
  //ztit=ztit(0, ztit.First('('));
  h2->SetTitle(Form("%s %.3f - %.3f;%s;%s",ztit.Data(), hh->GetZaxis()->GetBinLowEdge(iz0),hh->GetZaxis()->GetBinUpEdge(iz1), hh->GetXaxis()->GetTitle(), hh->GetYaxis()->GetTitle()));

  Double_t ntot = 0;
  for(Int_t ix=1; ix<= hh->GetNbinsX(); ix++){
    for(Int_t iy=1; iy<= hh->GetNbinsY(); iy++){
      const Int_t nn = hh->Integral(ix,ix,iy,iy, iz0,iz1);
      h2->SetBinContent(ix,iy,nn);

      ntot += nn;
    }
  }

  h2->SetEntries(ntot);
  return h2;
}

TChain * NeutrinoTools::InputFiles(const TString file, const TString tr, Char_t *dir)
{
  TChain *ch=new TChain(tr);

  if(file.Contains(".root"))
    ch->Add(file);
  else{
    ifstream fin(file);
    if(!fin){
      printf("NeutrinoTools::InputFiles file not found \n%s\n\n",file.Data()); exit(1);
    }

    TString buff;
    while(fin.good()){
      fin>>buff;
      if(buff!=""){
        if(dir){
          buff.Prepend(dir);
        }
        ch->Add(buff);
      }
    }
  }

  //const Int_t ent = ch->GetEntries(); //takes infinity time!!
  printf("\t%d trees!\n",ch->GetNtrees());

  return ch;
}

Double_t * NeutrinoTools::GetAxisArray(TAxis * aa)
{
  const Int_t nbin=aa->GetNbins();
  Double_t *bins = new Double_t[nbin+1];

  for(Int_t ii=0; ii<=nbin; ii++){
    bins[ii] = aa->GetBinUpEdge(ii);
  }

  return bins;
}

void NeutrinoTools::FitSlicesY(const TH2D *hh,
           TH1D *&hnor,
           TH1D *&hmpv,
           TH1D *&hwid,
           TH1D *&hres,
           TH1D *&hchi,
           const TString formula,
           const Double_t thres,
           TList *ll)
{
  const Int_t x0 = hh->GetXaxis()->GetFirst();
  const Int_t x1 = hh->GetXaxis()->GetLast();
  const Int_t y0 = hh->GetYaxis()->GetFirst();
  const Int_t y1 = hh->GetYaxis()->GetLast();

  const Int_t nx = hh->GetNbinsX();
  const Int_t ny = hh->GetNbinsY();
  const Double_t xmin = hh->GetXaxis()->GetXmin();
  const Double_t xmax = hh->GetXaxis()->GetXmax();
  const Double_t ymin = hh->GetYaxis()->GetXmin();
  const Double_t ymax = hh->GetYaxis()->GetXmax();

  hnor = new TH1D(Form("%s_%samp", hh->GetName(), formula.Data()), "", nx, xmin, xmax); if(ll){ll->Add(hnor);}
  hmpv = new TH1D(Form("%s_%smpv",hh->GetName(), formula.Data()), "", nx, xmin, xmax); if(ll){ll->Add(hmpv);}
  hwid = new TH1D(Form("%s_%swid",hh->GetName(), formula.Data()), "", nx, xmin, xmax); if(ll){ll->Add(hwid);}
  hres = new TH1D(Form("%s_%sres",hh->GetName(), formula.Data()), "", nx, xmin, xmax); if(ll){ll->Add(hres);}
  hchi = new TH1D(Form("%s_%schi",hh->GetName(), formula.Data()), "", nx, xmin, xmax); if(ll){ll->Add(hchi);}

  const Double_t *hxbins = GetAxisArray(hh->GetXaxis());
  hnor->GetXaxis()->Set(nx, hxbins);
  hmpv->GetXaxis()->Set(nx, hxbins);
  hwid->GetXaxis()->Set(nx, hxbins);
  hres->GetXaxis()->Set(nx, hxbins);
  hchi->GetXaxis()->Set(nx, hxbins);
  delete hxbins;

  for(Int_t ix=x0; ix<=x1; ix++){
    TH1D *htmp = new TH1D(Form("%s_%s%d", hh->GetName(), formula.Data(), ix),"",ny, ymin, ymax);
    //checked, ok
    const Double_t *hhybins = GetAxisArray(hh->GetYaxis());
    htmp->GetXaxis()->Set(ny, hhybins);
    delete hhybins;

    Double_t ntot = 0;
    for(Int_t iy=y0; iy<=y1; iy++){
      const Double_t be = hh->GetBinError(ix,iy);
      const Double_t bc = hh->GetBinContent(ix, iy);

      if(be<EPSILON){
        if(bc>EPSILON){
          printf("NeutrinoTools::FitSlicesY error %d %d %e %e\n", ix, iy, be, bc); exit(1);
        }
        continue;
      }

      htmp->SetBinContent(iy, bc);
      htmp->SetBinError(iy, be);

      ntot += (bc/be)*(bc/be);

      //if(be) printf("test %d %d : %f %f %f\n", ix, iy, bc, be, pow(bc/be,2));
    }


    hnor->SetBinContent(ix, ntot);
    hnor->SetBinError(  ix, TMath::Sqrt(ntot));

    if(ntot<thres || htmp->GetRMS()<EPSILON){
      delete htmp;
      continue;
    }

    //test htmp->Draw();
    Double_t pars[10]={htmp->Integral(0,htmp->GetNbinsX()+1)*htmp->GetBinWidth(1), htmp->GetMean()*0.9, htmp->GetRMS()*0.5};
    Double_t errs[10]={0,0,0,0,0,0,0,0,0,0}, chi[10]={0,0,0,0,0,0,0,0,0,0};

    if(formula.Contains("RMS")){
      pars[1]=htmp->GetMean();
      errs[1]=htmp->GetMeanError();

      //remember that GetRMS is affected by SetRangeUser!!
      pars[2]=htmp->GetRMS();
      errs[2]=htmp->GetRMSError();
    }
    else{
      TF1 * tmpf1 = 0x0;
      if(formula.Contains("Gaus")){
        tmpf1 = new TF1("tmpf1","TMath::Abs([0])*TMath::Gaus(x,[1],TMath::Abs([2]),1)", xmin, xmax);
      }
      else if(formula.Contains("Cauchy")){
        tmpf1 = new TF1("tmpf1","TMath::Abs([0])*TMath::CauchyDist(x,[1],TMath::Abs([2]))", xmin, xmax);
      }
      else{
        printf("NeutrinoTools::FitSlicesY not known formula %s\n", formula.Data());exit(1);
      }

      tmpf1->SetParameters(pars);
      htmp->Fit(tmpf1, Form("%sQ",formula.Contains("LS")?"":"L"));

      tmpf1->GetParameters(pars);
      for(Int_t ipar=0; ipar<tmpf1->GetNpar(); ipar++){
        errs[ipar]=tmpf1->GetParError(ipar);
      }

      chi[0]=tmpf1->GetChisquare();
      chi[1]=tmpf1->GetNDF();

      delete tmpf1;
    }

    pars[2]=TMath::Abs(pars[2]);
    //hnor->SetBinContent(ix, htmp->GetBinContent(htmp->GetMaximumBin()));//htmp->Integral(0,htmp->GetNbinsX()+1));
    hmpv->SetBinContent(ix, pars[1]);
    if(formula.Contains("+")){
      hmpv->SetBinError(  ix, pars[2]);
    }
    else{
      hmpv->SetBinError(  ix, errs[1]);
    }

    hwid->SetBinContent(ix, pars[2]);
    hwid->SetBinError(  ix, errs[2]);

    hres->SetBinContent(ix, fabs(pars[1])>EPSILON? pars[2]/fabs(pars[1]):0);
    hres->SetBinError(  ix, fabs(pars[1])>EPSILON? errs[2]/fabs(pars[1]):0);

    hchi->SetBinContent(ix, chi[1]>=1 ? chi[0]/chi[1]: 0);
    hchi->SetBinError(ix, 0);

    if(ll){
      ll->Add(htmp);
    }
    else{
      delete htmp;
    }
  }

  TH1 *hhs[]={hnor, hmpv, hwid, hres, hchi};
  const TString yt[]={"N", "Mean", "#sigma", "#sigma/MPV", "#chi^{2}/NDOF"};
  const Int_t nh = sizeof(hhs)/sizeof(TH1*);
  for(Int_t ii=0; ii<nh; ii++){
    hhs[ii]->SetYTitle(Form("%s of %s", yt[ii].Data(), hh->GetYaxis()->GetTitle()));
    hhs[ii]->SetXTitle(hh->GetXaxis()->GetTitle());
    hhs[ii]->GetYaxis()->SetTitleOffset(hh->GetYaxis()->GetTitleOffset());
    hhs[ii]->SetTitle(hh->GetTitle());
  }
}

void NeutrinoTools::ScaleToRef(TH1D * hh, const TH1D *href)
{
  for(Int_t ii=0; ii<=hh->GetNbinsX()+1; ii++){
    const Double_t den = href->GetBinContent(ii);

    Double_t val = 0;
    Double_t err = 0;
    if(den>EPSILON){
      val = hh->GetBinContent(ii)/den;
      err = hh->GetBinError(ii)/den;
    }

    hh->SetBinContent(ii, val);
    hh->SetBinError(ii,err);
  }
}

TGraphAsymmErrors * NeutrinoTools::GetFluxMap(const Double_t tol, const TH2D *hh, TH2D * & hmap, TH2D * & hdiff)
{
  //
  //tol: tolerance
  //use variable y to estimate x
  //chi2 not useful due to very small stat. err
  //

  if(hmap) delete hmap;
  hmap = (TH2D*) hh->Clone(Form("%smap", hh->GetName()));
  hmap->Scale(0);
  hmap->SetTitle(Form("%s Flux Map", hh->GetTitle()));

  if(hdiff) delete hdiff;
  hdiff= (TH2D*) hmap->Clone(Form("%sdiff", hmap->GetName()));
  hdiff->SetTitle(Form("%s Abs. Diff. to 1", hmap->GetTitle()));

  const Int_t nx = hmap->GetNbinsX();
  const Int_t ny = hmap->GetNbinsY();

  for(Int_t ix=1; ix<=nx; ix++){
    for(Int_t iy=1; iy<=ny; iy++){
      const Double_t Xintegral = hh->Integral(ix, nx+1, 0,  ny+1);
      if(Xintegral<EPSILON){
        continue;
      }

      const Double_t Yintegral = hh->Integral(0,  nx+1, iy, ny+1);

      const Double_t fluxratio = Yintegral/Xintegral;
      const Double_t frerr = TMath::Sqrt(fluxratio*(1-fluxratio)/Xintegral);

      hmap->SetBinContent(ix, iy, fluxratio);
      hmap->SetBinError(ix, iy, frerr);

      const Double_t frdiff = TMath::Abs(fluxratio-1);
      hdiff->SetBinContent(ix, iy, frdiff);
      hdiff->SetBinError(ix, iy, frerr);
    }
  }

  TGraphAsymmErrors * gr = new TGraphAsymmErrors;

  //use y to estimate x, so grX:= mapY
  for(Int_t iy=1; iy<=ny; iy++){
    Int_t varid = kIniValue;
    Double_t mindiff = 1e10;

    for(Int_t ix=1; ix<=nx; ix++){
      //xintegral<EPSILON, bin not set
      if(hdiff->GetBinError(ix,iy)<EPSILON){
        continue;
      }

      const Double_t frd = hdiff->GetBinContent(ix, iy);
      if( frd<mindiff ){
        varid = ix;
        mindiff = frd;
      }
    }

    if(mindiff<tol){
      const Double_t varFromCe = hdiff->GetYaxis()->GetBinCenter(iy);
      const Double_t varFromEl = varFromCe - hdiff->GetYaxis()->GetBinLowEdge(iy);
      const Double_t varFromEh = hdiff->GetYaxis()->GetBinUpEdge(iy) - varFromCe;

      const Double_t varToCe = hdiff->GetXaxis()->GetBinCenter(varid);
      const Double_t varToEl = varToCe - hdiff->GetXaxis()->GetBinLowEdge(varid);
      const Double_t varToEh = hdiff->GetXaxis()->GetBinUpEdge(varid) - varToCe;

      const Int_t ip = gr->GetN();
      gr->SetPoint(ip, varFromCe, varToCe);
      gr->SetPointError(ip, varFromEl, varFromEh, varToEl, varToEh);
    }
  }

  gr->SetName(Form("%scurve", hh->GetName()));
  gr->SetTitle(Form("half integrated (to infinity) flux map (tolerance %.2f%%)", tol));
  gr->GetXaxis()->SetTitle(hdiff->GetYaxis()->GetTitle());
  gr->GetYaxis()->SetTitle(hdiff->GetXaxis()->GetTitle());

  return gr;
}

Double_t NeutrinoTools::SampleFermiMomentum()
{
  //http://geant4.cern.ch/G4UsersDocuments/UsersGuides/PhysicsReferenceManual/html/node131.html#SECTION051212000000000000000
  /*
 For light nuclei with $A < 17$ nucleon density is given by a harmonic oscillator shell model [3], e. g.
\begin{displaymath}
\rho(r_i) = (\pi R^2)^{-3/2}\exp{(-r_i^2/R^2)},
\end{displaymath}	(25.3)

where  $R^2 = 2/3<r^2> = 0.8133 A^{2/3}$ fm$^2$. To take into account nucleon repulsive core it is assumed that internucleon distance $d > 0.8$ fm;
The nucleus is assumed to be isotropic, i.e. we place each nucleon using a random direction and the previously determined radius  $r_i$.
The initial momenta of the nucleons $p_i$ are randomly choosen between $0$ and $p^{max}_F(r)$, where the maximal momenta of nucleons (in the local Thomas-Fermi approximation [4]) depends from the proton or neutron density $\rho $ according to
\begin{displaymath}
p^{max}_F(r) = \hbar c(3\pi^2 \rho(r))^{1/3}
\end{displaymath}	(25.4)

To obtain momentum components, it is assumed that nucleons are distributed isotropic in momentum space; i.e. the momentum direction is chosen at random.
   */

  //C12
  const Double_t A1 = 12;
  const Double_t R2 = 0.8133*TMath::Power(A1,2./3.);
  const Double_t R1 = TMath::Sqrt(0.8133)*TMath::Power(A1,1./3.);

  //output NeutrinoTools::SampleFermiMomentum A1 12.000000 , R2 4.262898 , R1 2.064679//x-checked with wolframalpha
  //printf("NeutrinoTools::SampleFermiMomentum A1 %f , R2 %f , R1 %f\n", A1, R2, R1);

  //rr distribution (1000000 entries) tested with fit to exp(-x*x/[1]) and [1]=4.25248e+00, = R2
  const Double_t s2 = 1.41421356; //TMath::Sqrt(2);
  const Double_t r1 = fRan->Gaus(0 ,R1/s2);

  //Mev fm http://pdg.lbl.gov/2014/reviews/rpp2014-rev-phys-constants.pdf
  const Double_t l2e = 197.327;

  //http://www.wolframalpha.com/input/?i=3%5E%281%2F3%29*pi%5E%281%2F6%29
  const Double_t pmax = l2e * 1.7454151 * (1./R1) * TMath::Exp(-r1*r1/3./R2) / 1e3 ; //in GeV

  /*
    #include "NeutrinoTools.h"
    TH1D *hh = new TH1D("hh","",500,0,0.2)
    for(int ii=0; ii<1000000; ii++){double rr= NeutrinoTools::SampleFermiMomentum(); hh->Fill(rr);}
    hh->Draw()
   */

  return fRan->Rndm()*pmax;
}

void NeutrinoTools::SampleUnityIsotropicVector(TVector3 * vec)
{
  const Double_t phi = fRan->Rndm()*TMath::TwoPi();
  const Double_t rr = fRan->Rndm();
  const Double_t theta = TMath::ACos(1-2*rr);

  vec->SetMagThetaPhi(1,theta,phi);
}


void NeutrinoTools::Deflection(TVector3 * vec, const Double_t dtheta, const Double_t dphi)
{
  /*
//test

.L dotest.C+
TH3D *hh = new TH3D("hh","",100,-5,5,100,-5,5,100,-5,5)
hh->Draw()
TVector3 aa(1,1,1)
const TVector3 a0=aa;
const Int_t nn=50;
hh->Fill(aa.X(),aa.Y(),aa.Z());
for(int ii=0; ii<(3*nn); ii++){ Deflection(&aa, (TMath::PiOver4()),(TMath::TwoPi())/nn*ii); hh->Fill(aa.X(),aa.Y(),aa.Z()); aa=a0;}
hh->SetMarkerStyle(20);
hh->SetMarkerSize(1);
hh->Draw();

  */
  const Double_t theta0 = vec->Theta();
  const Double_t phi0 = vec->Phi();

  vec->SetMagThetaPhi(vec->Mag(), dtheta, dphi);

  vec->RotateY(theta0);
  vec->RotateZ(phi0);
}

Double_t NeutrinoTools::RutherfordTheta(const Double_t regpar)
{
  const Double_t rr = fRan->Rndm();

  //2) Rutherfordian
  //eq 8.46 Bielajew A
  //a->inf isotropic
  //a->0 delta at 0

  /*
//test
#include "NeutrinoTools.h"
TH1D * hh=new TH1D("hh","",100,0,3.2)
for(int ii=0; ii<1000000; ii++){double aa = NeutrinoTools::RutherfordTheta(1e-1);hh->Fill(aa);}
hh->Draw()
TF1 * f1 = new TF1("f1","[0]*sin(x)/pow(1-cos(x)+1e-1,2)",0,4)
f1->SetParameter(0,1)
hh->Fit(f1,"L")

   */
  const Double_t theta = TMath::ACos( 1-2*regpar*(1-rr)/(regpar+2*rr) );

  return theta;
}

Double_t NeutrinoTools::RutherfordEnergy(const Double_t regpar)
{
  const Double_t rr = fRan->Rndm();

  //f~1/E^2 -> a/(E+a)^2
  //a->inf: flat
  //a->0: delta at 0
  /*
//test
#include "NeutrinoTools.h"
TH1D * hh=new TH1D("hh","",100,0.1,30)
for(int ii=0; ii<1000000; ii++){double aa = NeutrinoTools::RutherfordEnergy(1e-1);hh->Fill(aa);}
hh->Draw()
TF1 * f1 = new TF1("f1","[0]/pow(x+1e-1,2)",0.1,40)
f1->SetParameter(0, 1)
hh->Fit(f1,"L")

   */

  return rr*regpar/(1-rr);
}

void NeutrinoTools::RutherfordTransport(TLorentzVector * v0, const Int_t nstep, const Double_t regphi, const Double_t regtheta, const Double_t regde)
{
  const Double_t mass = v0->M();

  for(Int_t ii=0; ii<nstep; ii++){

    const Double_t de = RutherfordEnergy(regde);

    const Double_t eNew = (v0->E()-de);
    if(eNew<mass){
      continue;
    }

    const Double_t pNew = TMath::Sqrt( eNew*eNew - mass*mass );

    //--

    const Double_t dphi = regphi*fRan->Rndm()*TMath::TwoPi();
    const Double_t dtheta = RutherfordTheta(regtheta);

    TVector3 dir = v0->Vect();
    Deflection(&dir, dtheta, dphi);

    const Double_t thetaNew = dir.Theta();
    const Double_t phiNew = dir.Phi();

    TVector3 vNew;
    vNew.SetMagThetaPhi(pNew,thetaNew, phiNew);
    v0->SetVectM(vNew, mass);
  }
}

Double_t NeutrinoTools::GetPtFast(const TVector3 * refdir, const Double_t dir0, const Double_t dir1, const Double_t dir2, const Double_t pp)
{
  TVector3 mom(dir0,dir1,dir2);

  if(mom.Mag()<EPSILON){
    printf("NeutrinoTools::GetPtFast mom.Mag null!\n");
    mom.Print();
    printf("\n");
    return kZeroDivider;
  }

  mom *= pp*1e-3/mom.Mag();

  return mom.Mag()*GetSin(refdir, &mom, "NeutrinoTools::GetPtFast");
}

const TVector3 *NeutrinoTools::GetVecT(const TLorentzVector * refdir, const Double_t xx, const Double_t yy, const Double_t zz)
{
  //
  //w.r.t. beam direction
  //
  if(!refdir){
    printf("TVector3 *NeutrinoTools::GetVecT refdir null\n"); exit(1);
  }

  const TVector3 vec(xx,yy,zz);

  TVector3 vRotated(vec);
  vRotated.Rotate(TMath::Pi(), refdir->Vect());

  const TVector3 *vt = new TVector3( (vec - vRotated)*0.5 );

  return vt;
}


Double_t NeutrinoTools::CalcAlpha(const TLorentzVector *refdir, const TLorentzVector * mup4)
{
  //
  //0~ 180
  //
  return GetAngle(refdir, mup4, "NeutrinoTools::CalcAlpha");
}

Double_t NeutrinoTools::GetCos(const TVector3 * v1, const TVector3 * v2, const TString tag)
{
  if(!v1 || !v2){
    printf("NeutrinoTools::GetCos v1 or v2 null!! %s\n", tag.Data()); exit(1);
  }

  if(v1->Mag()<EPSILON || v2->Mag()<EPSILON){
    printf("NeutrinoTools::GetCos v1 or v2 null!! %s\n", tag.Data());
    v1->Print();
    v2->Print();
    printf("\n");
    return kZeroDivider;
  }

  return v1->Dot(*v2) /v1->Mag()/v2->Mag();
}

Double_t NeutrinoTools::GetAngle(const TVector3 * v1, const TVector3 * v2, const TString tag)
{
  const Double_t tmpcos = GetCos(v1, v2, tag+"GetAngleGetCos");
  if(tmpcos== kZeroDivider){
    return kZeroDivider;
  }

  return TMath::ACos( tmpcos ) * TMath::RadToDeg();
}

Double_t NeutrinoTools::GetSin(const TVector3 * v1, const TVector3 * v2, const TString tag)
{
  const Double_t tmpcos = GetCos(v1, v2, tag+"GetSinGetCos");
  if(tmpcos== kZeroDivider){//including kZeroDivider
    return kZeroDivider;
  }

  if( (1-tmpcos*tmpcos) < EPSILON){
    printf("NeutrinoTools::GetSin tmpcos too close to 1 %f %f\n", tmpcos, 1-tmpcos*tmpcos); exit(1);
  }

  return TMath::Sqrt(1-tmpcos*tmpcos);

  /*
  if(!v1 || !v2){
    printf("NeutrinoTools::GetSin v1 or v2 null!!"); exit(1);
  }
  const TVector3 tmp = v1->Cross(*v2);

  if(v1->Mag()<EPSILON || v2->Mag()<EPSILON){
    printf("NeutrinoTools::GetSin v1 or v2 null!!");
    v1->Print();
    v2->Print();
    exit(1);
  }

  return tmp.Mag()/v1->Mag()/v2->Mag();
  */
}

Double_t NeutrinoTools::GetMass(const TLorentzVector * par1, const TLorentzVector * par2)
{
  Double_t invmass = -999;

  if(par1 && par2){
    TLorentzVector sumvec = *par1 + *par2;
    invmass = sumvec.M();
  }

  return invmass;
}

Double_t NeutrinoTools::GetDeltaPTT(const TVector3 &p1, const TVector3 &p2, const TVector3 & nuDir, const TVector3 & refDir)
{
  const Double_t ptt1 = GetPTT(p1, nuDir, refDir);
  const Double_t ptt2 = GetPTT(p2, nuDir, refDir);

  //both var obtained from GetPTT are signed in the same coordinate
  //const Double_t sign = ptt1>0? 1 : -1;
  const Double_t sign = 1;;//symmetrized if only Hydrogen is concerned!

  return (ptt1+ptt2)*sign;
}

Double_t NeutrinoTools::GetPTT(const TVector3 & parMom, const TVector3 & nuDir, const TVector3 & muonDir, const Bool_t kprint)
{
  const TVector3 newaxis = (nuDir.Cross(muonDir)).Unit();

  const Double_t ptt = newaxis.Dot(parMom);

  if(kprint){
    parMom.Print();
    nuDir.Print();
    muonDir.Print();
    newaxis.Print();
    printf("NeutrinoTools::GetPTT ptt %f\n", ptt);
  }

  return ptt;
}

Double_t NeutrinoTools::GetCosDeltaPhiT(const TVector3 * pt1, const TVector3 * pt2)
{
  return pt1->Dot(*pt2)*(-1)/pt1->Mag()/pt2->Mag();
}

/*
//the directionality w.r.t. the plane is lost, use PTT instead
Double_t NeutrinoTools::GetKT(const TVector3 * ptScal, const TVector3 * ptAxis)
{
  const Double_t arg = 1-TMath::Power(GetCosDeltaPhiT(ptScal, ptAxis),2);
  if(arg<EPSILON){
    return 0;
  }
  else{
    return ptScal->Mag()*TMath::Sqrt(arg);
  }
}
*/

void NeutrinoTools::SetDeltaPt(TVector3 * deltapt, const TVector3 * ptmuon, const TVector3 * ptproton)
{
  //ptmuon and ptproton already in the same plain which is perpendicular to the neutrino and already in a near back-to-back configuration
  const TVector3 tmpd = (*ptmuon)+(*ptproton);

  const Double_t phi = TMath::ACos( GetCosDeltaPhiT(ptmuon, ptproton) );

  const Double_t theta = TMath::ACos( tmpd.Dot(*ptmuon)*(-1)/tmpd.Mag()/ptmuon->Mag()  );

  deltapt->SetMagThetaPhi(tmpd.Mag(),theta, phi);
}

Double_t NeutrinoTools::GetEnuApp(const TVector3 & v3Nu, const Double_t mN, const TLorentzVector *lvL, const Double_t mX)
{
  const Double_t eL = lvL->E();
  const TVector3 pL = lvL->Vect();

  const TVector3 v3U = v3Nu.Unit();

  const Double_t deno = 2*( eL-mN-pL.Dot(v3U) )+EPSILON;
  const Double_t enu = ( TMath::Power(mN-eL, 2) -pL*pL - mX*mX)/deno;

  return enu;

  /*
//test
#include "NeutrinoTools.h"
.L AnaTask1.C+
const TVector3 tmpnu(0,0,1);
const Double_t tmpn=NeutrinoTools::NeutronMass();
const Double_t tmpx=NeutrinoTools::ProtonMass();
const TLorentzVector tmpmu(0,0,0,NeutrinoTools::MuonMass());
Double_t aa = NeutrinoTools::GetEnuApp(tmpnu, tmpn, tmpmu, tmpx);
//root [7] aa
//(Double_t)1.10896155380221545e-01

#include "NeutrinoTools.h"
.L AnaTask1.C+
const TVector3 tmpnu(0,0,1);
const Double_t tmpn=NeutrinoTools::ProtonMass();
const Double_t tmpx=NeutrinoTools::DeltaPPMass();
const TLorentzVector tmpmu(0,0,0,NeutrinoTools::MuonMass());
Double_t aa = NeutrinoTools::GetEnuApp(tmpnu, tmpn, &tmpmu, tmpx)
aa
//(Double_t)4.95174780611457477e-01

//electron
const TLorentzVector tmpmu(0,0,0,0.511/1e3);
Double_t aa = NeutrinoTools::GetEnuApp(tmpnu, tmpn, &tmpmu, tmpx)
aa
//(Double_t)3.40400329550568093e-01

#include "NeutrinoTools.h"
.L AnaTask1.C+
const TLorentzVector gnu(0,0,5,5);
const TVector3 tmpnu=gnu.Vect();

const TLorentzVector gN(0,0,0,NeutrinoTools::NeutronMass());
const Double_t tmpn=gN.M();

TLorentzVector tmpmu; tmpmu.SetXYZM(0.3,0.4,0.5,NeutrinoTools::MuonMass());

const Double_t tmpx=(gnu+gN-tmpmu).M()

Double_t aa = NeutrinoTools::GetEnuApp(tmpnu, tmpn, tmpmu, tmpx);
aa-gnu.E()
   */
}

Int_t NeutrinoTools::FillTH1I(TList *lin, const Int_t id, const Int_t var)
{
  //printf("test NeutrinoTools::FillTH %d %d\n", id, var);
  return ((TH1I*) lin->At(id))->Fill(var);
}

Int_t NeutrinoTools::FillTH1D(TList *lin, const Int_t id, const Double_t var)
{
  //printf("test NeutrinoTools::FillTH %d %f\n", id, var);
  return ((TH1D*) lin->At(id))->Fill(var);
}

Int_t NeutrinoTools::FillTH2D(TList *lin, const Int_t id, const Double_t var1, const Double_t var2)
{
  //printf("test NeutrinoTools::FillTH %d %f %f\n", id, var1, var2);
  return ((TH2D*) lin->At(id))->Fill(var1, var2);
}

Int_t NeutrinoTools::FillTH3D(TList *lin, const Int_t id, const Double_t var1, const Double_t var2, const Double_t var3)
{
  //printf("test NeutrinoTools::FillTH %d %f %f\n", id, var1, var2, var3);
  return ((TH3D*) lin->At(id))->Fill(var1, var2, var3);
}

TString NeutrinoTools::GetTitleFromVar(const TString var)
{
  TString tit(var);
  tit.ReplaceAll("fNeutrinoSim->E()","E_{#nu sim} (GeV)");
  tit.ReplaceAll("fQ2Sim","Q^{2}_{sim} (GeV^{2}/c^{4})");
  tit.ReplaceAll("fMultiplicity","N_{all rec}");
  tit.ReplaceAll("fNeutrinoParentType","1:e, 2:#mu, 3:K_{L}, 4:#pi^{#pm}, 5:K^{#pm}, 6:p");
  tit.ReplaceAll("fNeutrinoRec->Theta()*TMath::RadToDeg()","#theta_{#nu rec} (deg)");
  tit.ReplaceAll("fNeutrinoSim->Theta()*TMath::RadToDeg()","#theta_{#nu sim} (deg)");
  tit.ReplaceAll("(fNeutMode==1||fNeutMode==2)","CCQE");

  tit.ReplaceAll("fDeltaRecPt->Mag()","#delta p_{T}^{rec} (GeV/c)");
  tit.ReplaceAll("fDeltaSimPt->Mag()","#delta p_{T}^{sim} (GeV/c)");

  tit.ReplaceAll("fDeltaRecPt->Phi()","#delta#phi_{T}^{rec} (rad)");
  tit.ReplaceAll("fDeltaSimPt->Phi()","#delta#phi_{T}^{sim} (rad)");

  tit.ReplaceAll("fDeltaRecPt->Theta()","#delta#alpha_{T}^{rec} (rad)");
  tit.ReplaceAll("fDeltaSimPt->Theta()","#delta#alpha_{T}^{sim} (rad)");

  tit.ReplaceAll("fProtonPTTRec","k_{T}^{p rec} (GeV/c)");
  tit.ReplaceAll("fProtonPTTSim","k_{T}^{p sim} (GeV/c)");

  tit.ReplaceAll("fMuonKTRec", "k_{T}^{#mu rec} (GeV/c)");
  tit.ReplaceAll("fMuonKTSim", "k_{T}^{#mu sim} (GeV/c)");

  //=======================

  tit.ReplaceAll("fMomErr/fMuonRec->P()","#sigma_{#mu p}/p rec");
  tit.ReplaceAll("fMuonSim->P()","p_{#mu sim} (GeV/c)");
  tit.ReplaceAll("fMuonSim->Pt()","p_{t}^{#mu sim int} (GeV/c)");
  tit.ReplaceAll("fMuonRec->Pt()","p_{t}^{#mu rec int} (GeV/c)");
  //tit.ReplaceAll("fMuonPtGuess->Mag()","p_{t}^{#mu rec, #nu sim} (GeV/c)");
  tit.ReplaceAll("fMuonSimPt->Mag()","p_{t}^{#mu sim} (GeV/c)");
  tit.ReplaceAll("fMuonRec->P()","p_{#mu rec} (GeV/c)");
  tit.ReplaceAll("fMuonRecPt->Mag()","p_{t}^{#mu rec} (GeV/c)");
  tit.ReplaceAll("fMuonRecFlightPath->Pt()","#mu lever arm (m)");
  tit.ReplaceAll("fMuonRecFlightPath->Mag()","#mu straight path length (m)");
  tit.ReplaceAll("(fMuonCharge[1]==-1)*2+(fMuonCharge[0]==-1)","#mu TPC charge (11,-11,1-1,-1-1)");
  tit.ReplaceAll("fMuonRecNhits","selmu_tpc_nhits");
  tit.ReplaceAll("fMuonTypeSim","#mu ID (1:e, 2:#mu, 3:K_{L}, 4:#pi^{#pm}, 5:K^{#pm}, 6:p)");
  tit.ReplaceAll("fMuonAlphaRec","#theta_{#mu rec} (deg)");
  tit.ReplaceAll("fMuonAlphaSim","#theta_{#mu sim} (deg)");
  tit.ReplaceAll("fMuonRecEndPos->X()","#mu endposX (m)");
  tit.ReplaceAll("fMuonRecEndPos->Y()","#mu endposY (m)");
  tit.ReplaceAll("fMuonRecEndPos->Z()","#mu endposZ (m)");
  tit.ReplaceAll("fMuonRecVertex->X()","#mu vertexX (m)");
  tit.ReplaceAll("fMuonRecVertex->Y()","#mu vertexY (m)");
  tit.ReplaceAll("fMuonRecVertex->Z()","#mu vertexZ (m)");
  tit.ReplaceAll("fMuonRec->Theta()*TMath::RadToDeg()","#theta_{#mu rec lab} (deg)");
  tit.ReplaceAll("fMuonRecChi2/fMuonRecNDOF","#chi^{2}_{#mu track}/NDOF");

  tit.ReplaceAll("fMomErr/fProtonRec->P()","#sigma_{p p}/p rec");
  tit.ReplaceAll("fProtonSim->P()","p_{p sim} (GeV/c)");
  tit.ReplaceAll("fProtonSim->Pt()","p_{t}^{p sim int} (GeV/c)");
  tit.ReplaceAll("fProtonRec->Pt()","p_{t}^{p rec int} (GeV/c)");
  //tit.ReplaceAll("fProtonPtGuess->Mag()","p_{t}^{p rec, #nu sim} (GeV/c)");
  tit.ReplaceAll("fProtonSimPt->Mag()","p_{t}^{p sim} (GeV/c)");
  tit.ReplaceAll("fProtonRec->P()","p_{p rec} (GeV/c)");
  tit.ReplaceAll("fProtonRecPt->Mag()","p_{t}^{p rec} (GeV/c)");
  tit.ReplaceAll("fProtonRecFlightPath->Pt()","p lever arm (m)");
  tit.ReplaceAll("fProtonRecFlightPath->Mag()","p straight path length (m)");
  tit.ReplaceAll("(fProtonCharge[1]==-1)*2+(fProtonCharge[0]==-1)","p TPC charge (11,-11,1-1,-1-1)");
  tit.ReplaceAll("fProtonRecNhits","selp_tpc_nhits");
  tit.ReplaceAll("fProtonTypeSim","p ID (1:e, 2:#mu, 3:K_{L}, 4:#pi^{#pm}, 5:K^{#pm}, 6:p)");
  tit.ReplaceAll("fProtonAlphaRec","#theta_{p rec} (deg)");
  tit.ReplaceAll("fProtonAlphaSim","#theta_{p sim} (deg)");
  tit.ReplaceAll("fProtonRecEndPos->X()","p endposX (m)");
  tit.ReplaceAll("fProtonRecEndPos->Y()","p endposY (m)");
  tit.ReplaceAll("fProtonRecEndPos->Z()","p endposZ (m)");
  tit.ReplaceAll("fProtonRecVertex->X()","p vertexX (m)");
  tit.ReplaceAll("fProtonRecVertex->Y()","p vertexY (m)");
  tit.ReplaceAll("fProtonRecVertex->Z()","p vertexZ (m)");
  tit.ReplaceAll("fProtonRec->Theta()*TMath::RadToDeg()","#theta_{p rec lab} (deg)");
  tit.ReplaceAll("fProtonRecChi2/fProtonRecNDOF","#chi^{2}_{p track}/NDOF");

  return tit;
}

/*
TString NeutrinoTools::GetPIDCut()
{
  //return "( f@RecNhits[0]<80 && f@RecFlightPath->Mag()>2.5 && f@Charge[0]==-1 && f@Charge[1]==-1 && fMultiplicity<4)";
  return "( f@RecNhits[0]<80 && f@RecFlightPath->Mag()>3.0 )";
}

TString NeutrinoTools::GetPtResCut()
{
  //don't cut on p, at least not on lower end, it is physics; too large P is mostly fluctuation
  //return "( f@Rec->P()<3.5 && f@RecNhits[0]>60 && f@RecFlightPath->Pt()>1.5  )";
  return "( 0.25<(f@RecChi2/f@RecNDOF) && (f@RecChi2/f@RecNDOF)<4 && f@RecNhits[0]>60 && f@RecFlightPath->Pt()>1.5 )";
}


TString NeutrinoTools::GetQACut()
{
  return Form("( %s && %s )", GetPIDCut().Data(), GetPtResCut().Data());
}

*/

TString NeutrinoTools::GetModeCut(const Int_t imode, TString & name)
{
  const TString modecuts[]={"1",
                            "fNeutrinoType==14 && ( fNeutMode==1 || fNeutMode==2 )",
                            "fNeutrinoType==14 && ( fNeutMode==11 || fNeutMode==12 || fNeutMode==13 || fNeutMode== 16 || fNeutMode==17 || fNeutMode== 22 || fNeutMode== 23 )",
                            "fNeutrinoType==14 && ( fNeutMode==21 )",
                            "fNeutrinoType==14 && ( fNeutMode==26 )",
                            //"fNeutrinoType!=14 && ( abs(fNeutMode)<30 )"
                            "(fNeutrinoType!=14 || abs(fNeutMode)>30)"
  };
  const TString modenames[]={"incl","#nu_{#mu}CCQE","#nu_{#mu}CC#Delta","#nu_{#mu}CCm#pi","#nu_{#mu}CCDIS","#slash{#nu_{#mu}} or NC"};

  name=modenames[imode];
  return Form("(%s)",modecuts[imode].Data());
}
