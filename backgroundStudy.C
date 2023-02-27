#include <stdio.h>
#include <TROOT.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH2F.h>
#include <Riostream.h>
#include <TRandom.h>
#include <TVector3.h>
#include <TVector2.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TF1.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TStyle.h>


// get Ring-radius from Cherenkov-angle
double getRadiusFromCkov(double ckovAngle);
Float_t GetFreonIndexOfRefraction(Float_t x);
Float_t GetQuartzIndexOfRefraction(Float_t x);
Double_t BackgroundFunc(Double_t *x, Double_t *par);

const double defaultPhotonEnergy = 6.75; 
const double refIndexFreon = GetFreonIndexOfRefraction(defaultPhotonEnergy);
const double refIndexQuartz = GetQuartzIndexOfRefraction(defaultPhotonEnergy);
const double  refIndexCH4 = 1.00; 




const double CH4GapWidth = 8;
const double  RadiatorWidth = 1.;
const double  QuartzWindowWidth = 0.5;
const double  EmissionLenght = RadiatorWidth/2;
  

void backgroundStudy(Int_t NumberOfEvents, Int_t NumberOfClusters, Double_t Hwidth)   
{

  const auto numberOfCkovPhotons = gRandom->Poisson(12);

  Double_t ThetaP=0,PhiP=0,PhiF=0,DegThetaP=0,DegPhiP=0;
  //Float_t /*RadiatorWidth,*/ QuartzWindowWidth,CH4GapWidth,EmissionLenght;
  Float_t FreonIndexOfRefraction,QuartzIndexOfRefraction,CH4IndexOfRefraction;
  Double_t ThetaF1,ThetaF2,ThetaF=0,ThetaLimite;
  Float_t Xpi=0,Ypi=0,Xf=0,Yf=0,Xf1=0,Yf1=0,Xf2=0,Yf2=0,Xp=0,Yp=0; 

  Float_t PhotonEnergy = 6.75; 
  
  FreonIndexOfRefraction = GetFreonIndexOfRefraction(PhotonEnergy);
  QuartzIndexOfRefraction = GetQuartzIndexOfRefraction(PhotonEnergy);
  CH4IndexOfRefraction = 1.00;
  

  TH1F *ThetaHough = new TH1F("ThetaHough","ThetaHough",1000,0,1);
  TH1F *NPhoton    = new TH1F("NPhoton","NPhoton",20,0.,20.); 
  TH1F *NHough     = new TH1F("NHuogh","NHough",1000,0.,1.); 

  TH1F *hTheta     = new TH1F("Background","Background; angle [rad]; counts/1 mrad",750,0.,0.75);


  TH1F *hTheta2    = new TH1F("Photon Cherenkov Angle Distribution","Photon Cherenkov Angle Distribution;angle [rad]; counts/1 mrad",750,0.,0.75);
  TH1F *hThetawg   = new TH1F("Track Cherenkov Angle Distribution","Track Cherenkov Angle Distribution;angle [rad]; counts/1 mrad",750,0.,0.75);
  TH1F *hThetaRing = new TH1F("hThetaRing","Ring Cherenkov; Cherenkov Angle [rad];",1000,0.,1.);

  TH2F *hClusterMap = new TH2F("Cluster Map", "Cluster Map; x [Pads]; y [Pads]",1000,-10.,10.,1000,-10.,10.);
  TH2F *hSignalAndNoiseMap = new TH2F("Cluster Map2", "Cluster Map2; x [Pads]; y [Pads]",1000,-25.,25.,1000,-25.,25.);

  TH2F *hphotonMap = new TH2F("Photon Map","Photon Map; x [Pads]; y [Pads]",1000,-10.,10.,1000,-10.,10.);




  TH2F *hphotonMap2 = new TH2F("Photon Map2","Photon Map2; x [Pads]; y [Pads]",1000,-25.,25.,1000,-25.,25.);
  TH2F *signalMap = new TH2F("signalHitMap","signalHitMap",1000,0.,-30.,1000,0.,30.);
  TH2F *noiseMap = new TH2F("noiseHitMap","noiseHitMap",1000,0.,-30.,1000,0.,30.);


   
  Float_t Deltax = (RadiatorWidth+QuartzWindowWidth+CH4GapWidth-EmissionLenght)*TMath::Tan(ThetaP)*TMath::Cos(PhiP);
  Float_t Deltay = (RadiatorWidth+QuartzWindowWidth+CH4GapWidth-EmissionLenght)*TMath::Tan(ThetaP)*TMath::Sin(PhiP);
	  
  Xpi = Xp - Deltax;
  Ypi = Yp - Deltay;
	  
  Float_t ThetaCherenkov[100000] = {0x0}, PhiCherenkov[100000] = {0x0}, DegPhiCherenkov[100000] = {0x0};
    
  for(Int_t iEvt = 0; iEvt<NumberOfEvents; iEvt++){
    
    Printf("event number = %i",iEvt);
    
    Float_t Xcen[100000],Ycen[100000];
     
    DegThetaP = 4.;//0.*(1 - 2*gRandom->Rndm(iEvt));
    DegPhiP   = 360*gRandom->Rndm(iEvt);
        
    ThetaP = TMath::Pi()*DegThetaP/180;
    PhiP = TMath::Pi()*DegPhiP/180;  
     
    for(Int_t n1=0;n1<NumberOfClusters; n1++) {// clusters loop
      
    //  Printf("cluster = %i",n1);
      
      Xcen[n1] = 60*(1 - 2*gRandom->Rndm(n1));
      Ycen[n1] = 60*(1 - 2*gRandom->Rndm(n1));
      hphotonMap->Fill(Xcen[n1], Ycen[n1]);
      hphotonMap2->Fill(Xcen[n1], Ycen[n1]);



      noiseMap->Fill(Xcen[n1], Ycen[n1]);
      hSignalAndNoiseMap->Fill(Xcen[n1], Ycen[n1]);
            

//      Xcen[n1] = 130*gRandom->Rndm(n1); Ycen[n1] = 130*gRandom->Rndm(n1); 
      
      TVector3 v2(Xcen[n1]-Xpi-EmissionLenght*TMath::Tan(ThetaP)*TMath::Cos(PhiP),Ycen[n1]-Ypi-EmissionLenght*TMath::Tan(ThetaP)*TMath::Sin(PhiP),RadiatorWidth+QuartzWindowWidth+CH4GapWidth-EmissionLenght); 
      
      PhiF = v2.Phi();     



     // R_ckov = sin (theta_ckov)*(RadiatorWidth - EmissionLenght)
     // R_qz = sin(theta_qz) * QuartzWindowWidth
     // R_0 = sin(theta_0)*CH4GapWidth
     
	




     //theta_qz = TMath::Asin(sin_ckov*n_GAP/n_quartz);

     // sin_ckov * n_Radiator = sin_qz * n_quartz
     // sin_qz * n_quartz = sin_theta0 * n_GAP
       
      
      ThetaLimite = TMath::ASin(CH4IndexOfRefraction/QuartzIndexOfRefraction);
      
      Double_t ThetaF0 = TMath::ASin(QuartzIndexOfRefraction/FreonIndexOfRefraction*TMath::Sin(ThetaLimite))-0.00001;


      
    //  Printf("ThetaF0 = %f",ThetaF0*TMath::RadToDeg());
	      
      Double_t ThetaF01 = TMath::ASin((FreonIndexOfRefraction/QuartzIndexOfRefraction)*(TMath::Sin(ThetaF0)));      
      
      Double_t ThetaF02 = TMath::ASin((QuartzIndexOfRefraction/CH4IndexOfRefraction)*(TMath::Sin(ThetaF01)));
	      
      Float_t X01 = EmissionLenght*TMath::Tan(ThetaP)*TMath::Cos(PhiP);
	      
      Float_t Y01 =  EmissionLenght*TMath::Tan(ThetaP)*TMath::Sin(PhiP);
	      
      Float_t X02 = (RadiatorWidth - EmissionLenght)*TMath::Tan(ThetaF0)*TMath::Cos(PhiF)+QuartzWindowWidth*TMath::Tan(ThetaF01)*TMath::Cos(PhiF)+CH4GapWidth*TMath::Tan(ThetaF02)*TMath::Cos(PhiF);
	      
      Float_t Y02 = (RadiatorWidth - EmissionLenght)*TMath::Tan(ThetaF0)*TMath::Sin(PhiF) + QuartzWindowWidth*TMath::Tan(ThetaF01)*TMath::Sin(PhiF) + CH4GapWidth*TMath::Tan(ThetaF02)*TMath::Sin(PhiF);  
	  
      Float_t X0 = X01 + X02;
      Float_t Y0 = Y01 + Y02;
	      
      Double_t ThetaMin = 0;
//      Double_t ThetaMax = 0.75+ThetaP;
      Double_t ThetaMax = ThetaF0;
	      
      Xf = 999;
      Yf = 999;
      
      Int_t nWhile = 0;
      
      while(TMath::Sqrt((Xf-Xcen[n1]+Xpi)*(Xf-Xcen[n1]+Xpi)+(Yf-Ycen[n1]+Ypi)*(Yf-Ycen[n1]+Ypi))>0.0001)
		
	{ 
          nWhile++;
          
	  ThetaF = (Double_t) (0.5*(ThetaMax - ThetaMin) + ThetaMin);
	  
	  ThetaF1 = TMath::ASin((FreonIndexOfRefraction/QuartzIndexOfRefraction)*(TMath::Sin(ThetaF)));     
	  
	  ThetaF2 = TMath::ASin((QuartzIndexOfRefraction/CH4IndexOfRefraction)*(TMath::Sin(ThetaF1)));
	  
	  Xf1 = EmissionLenght*TMath::Tan(ThetaP)*TMath::Cos(PhiP);
		  
	  Yf1 =  EmissionLenght*TMath::Tan(ThetaP)*TMath::Sin(PhiP);
	  
	  Xf2 = (RadiatorWidth - EmissionLenght)*TMath::Tan(ThetaF)*TMath::Cos(PhiF)+QuartzWindowWidth*TMath::Tan(ThetaF1)*TMath::Cos(PhiF)+CH4GapWidth*TMath::Tan(ThetaF2)*TMath::Cos(PhiF);
	  
		  
	  Yf2 = (RadiatorWidth - EmissionLenght)*TMath::Tan(ThetaF)*TMath::Sin(PhiF) + QuartzWindowWidth*TMath::Tan(ThetaF1)*TMath::Sin(PhiF) + CH4GapWidth*TMath::Tan(ThetaF2)*TMath::Sin(PhiF);  
		  

	 

	  Xf = Xf1 + Xf2;
	  Yf = Yf1 + Yf2;
	  Printf("Positions Xf1 %f Yf1 %f Xf2 %f Yf2 %f" , Xf1, Yf1, Xf2, Yf2);
		  


	  if(TMath::Sqrt((Xf-X0)*(Xf-X0)+(Yf-Y0)*(Yf-Y0))>TMath::Sqrt((Xcen[n1]-Xpi-X0)*(Xcen[n1]-Xpi-X0)+(Ycen[n1]-Ypi-Y0)*(Ycen[n1]-Ypi-Y0)))
		    
	    {
	      ThetaMin = ThetaF;
	    }
		  
	  else 
		    
	    {
	      ThetaMax = ThetaF;   
	    }  
            
            if(nWhile>30) break;
            
	} // while 
        
	      
      TVector3 vP((TMath::Sin(ThetaP))*(TMath::Cos(PhiP)),(TMath::Sin(ThetaP))*(TMath::Sin(PhiP)),(TMath::Cos(ThetaP)));
      TVector3 vz(0.,0.,1.);
	      
      TVector3 v1 = vP.Cross(vz);
	      
      TVector3 vF((TMath::Sin(ThetaF))*(TMath::Cos(PhiF)),(TMath::Sin(ThetaF))*(TMath::Sin(PhiF)),(TMath::Cos(ThetaF)));
      
      if(ThetaP==0)
	{	      
	  ThetaCherenkov[n1] = ThetaF;		  
	  PhiCherenkov[n1] = PhiF;	      
	}	  
      else		
	{
	  vF.Rotate(ThetaP,v1);
	  
	  ThetaCherenkov[n1] = vF.Theta();
	  PhiCherenkov[n1] = vF.Phi();
	}
      
      DegPhiCherenkov[n1] = 180*PhiCherenkov[n1]/(TMath::Pi());

      if(DegPhiCherenkov[n1]<0) DegPhiCherenkov[n1]+=360;	      

      TVector2 v5(Xcen[n1]-Xp,Ycen[n1]-Yp);
      
      //Double_t Phi = v5.Phi();
      
      //Float_t DegPhi = 180*Phi/TMath::Pi(); 
      
      hTheta->Fill(ThetaCherenkov[n1]);
	      	      
      //if(k1==2233) cout << " ThetaCherenkov " <<ThetaCherenkov[n1] <<"   PhiCherenkov = "<<DegPhiCherenkov[n1]<<"  event number = "<<k1<< endl;

      Int_t dstep = (Int_t)DegPhiCherenkov[n1]/20;
      if(dstep==18) dstep = 0;
      if(dstep<0 || dstep>18) continue;	      
      
   }//clusters loop 	  
  
  Int_t HCS[700] = {0x0};
	  
  for(Int_t k2=0;k2<700;k2++)
    {
      
      HCS[k2] = 0;             
      Int_t NphotHough = 0;
      
      for(Int_t p=0;p<NumberOfClusters;p++)
	{
		  
	  if(ThetaCherenkov[p]>(0.001*k2) && ThetaCherenkov[p]<(0.001*Hwidth+0.001*k2)) NphotHough++;
	}
	      
      HCS[k2] = NphotHough;
	      
      NHough->Fill(0.001*k2+0.001/2.,(Float_t)NphotHough);
	      
   }
	  
  Int_t LocPos = TMath::LocMax(700,HCS);
	  
  Int_t NphotTot = 0;
  Float_t MeanTot = 0;
  
  for(Int_t p=0;p<NumberOfClusters;p++)
    {
      if(ThetaCherenkov[p]>(0.001*LocPos) && ThetaCherenkov[p]<(Hwidth+0.001*LocPos)) 
	{
	  ThetaHough->Fill(ThetaCherenkov[p]);
		  
	  NphotTot++;
	  MeanTot+=ThetaCherenkov[p];		 
	  // 	  ThetaPhoton[NphotTot-1]=ThetaCherenkov[p];
	}
    }
	    
   Float_t RingThetaCherenkov = MeanTot/(Float_t)NphotTot;

   NPhoton->Fill(NphotTot);
   
   if(NphotTot<3) continue;
   
   hThetaRing->Fill(RingThetaCherenkov); // gjor denne på weighted m signal
        
 }// events loop

 // endre weight til fra Recon fra Giacomo sin branch

 TF1 *fBackGround = new TF1("fBackGround",BackgroundFunc,0,0.75,1);
 
 TCanvas *c1 = new TCanvas("c1","c1",800,800);
 TCanvas *c2 = new TCanvas("c2","c2",800,800);
 TCanvas *c3 = new TCanvas("c3","c3",800,800);


 TCanvas *photonCanvas = new TCanvas("Photon Hit Map","Photon Hit Map",800,800);
 TCanvas *photonCanvas2 = new TCanvas("Photon Hit Map2","Photon Hit Map2",800,800);

 TCanvas *clusterCanvas = new TCanvas("Cluster Hit Map","Cluster Hit Map",800,800);
 TCanvas *clusterCanvas2 = new TCanvas("Cluster Hit Map2","Cluster Hit Map2",800,800);

 photonCanvas->cd();
 hphotonMap->Draw();


 photonCanvas2->Divide(2,1);
 photonCanvas2->cd(1);
 hphotonMap->Draw();

 photonCanvas2->cd(2);
 hphotonMap2->Draw();

 clusterCanvas->cd();
 hClusterMap->Draw();


 clusterCanvas2->Divide(2,1);
 clusterCanvas2->cd(1);
 hClusterMap->Draw();

 clusterCanvas2->cd(2);
 hSignalAndNoiseMap->Draw();

 



//signalMap noiseMap hSignalAndNoiseMap


 TCanvas *c4 = new TCanvas("c4","c4",600,400);
 
 c1->cd();
 hTheta->Fit(fBackGround,"RQ");
 hTheta->Draw();
 
 Printf("par0 = %f",fBackGround->GetParameter(0));
 
 for(Int_t i=0; i < numberOfCkovPhotons; i++) {
   
   double ckovAngle = gRandom->Gaus(0.5,0.012);		    // random CkovAngle
   double ringRadius = getRadiusFromCkov(ckovAngle); // R in photon-map
						

   double alpha = static_cast<double>((3.14159)*(1-2*gRandom->Rndm(1)));    // angle in photon-map (-pi to pi)

   // get x and y values of Photon-candidate:
   double x = static_cast<double>(TMath::Cos(alpha)*ringRadius);
   double y = static_cast<double>(TMath::Sin(alpha)*ringRadius);  


   signalMap->Fill(x,y);//signalMap noiseMap hSignalAndNoiseMap
   hSignalAndNoiseMap->Fill(x,y);

   hTheta2->Fill(ckovAngle);
   hThetawg->Fill(ckovAngle);
  } 
 
 for(Int_t n2=0;n2<NumberOfClusters; n2++) {
  
   if(ThetaCherenkov[n2]>0.75) continue;
   
   Int_t bin = (Int_t)(ThetaCherenkov[n2]/(0.001))+1;
   
  // Printf("bin = %i, hTheta->GetBinContent(bin)=%f",bin, hTheta->GetBinContent(bin));
   
  // if(hTheta->GetBinContent(bin)==0) bin=bin+1;
   
   Double_t weight = 1 - fBackGround->Eval(ThetaCherenkov[n2])/hTheta->GetBinContent(bin);

   hTheta2->Fill(ThetaCherenkov[n2]);
    
   hThetawg->Fill(ThetaCherenkov[n2],weight);
 }
 

 TCanvas *signalCanvas = new TCanvas("SignalMap","SignalMap", 800, 800);
 signalCanvas->Divide(2,1);
 signalCanvas->cd(1);
 signalMap->Draw();

 TCanvas *noiseCanvas = new TCanvas("Noise Map","Noise Map", 800, 800);
 noiseCanvas->Divide(2,1);
 noiseCanvas->cd(1);
 noiseMap->Draw();

 TCanvas *signalNoiseCanvas = new TCanvas("Signal and Noise Map","Signal and Noise Map", 800, 800);
 signalNoiseCanvas->Divide(2,1);
 signalNoiseCanvas->cd(1);
 hSignalAndNoiseMap->Draw();



 c2->cd();
 hTheta2->Draw(); 
 
 c3->cd();
 hThetawg->Draw();
 
 c4->cd();
 NPhoton->Draw();

 TCanvas* ringCanvas = new TCanvas("Ring Cherenkov", "Ring Cherenkov", 800,800);
 ringCanvas->cd();
 hThetaRing->Draw();

  gStyle->SetOptStat("erm");
 {
   TCanvas* testCanvas = new TCanvas("test" ,"test", 800,800);
   testCanvas->Divide(2,1);
   auto pad = static_cast<TPad *>(testCanvas->cd(1));
   //pad->SetLogy(1);
   hThetawg->Draw();

   testCanvas->cd(2);
   auto pad2 = static_cast<TPad *>(testCanvas->cd(2));
   hphotonMap->Draw();

 }

 {
   TCanvas* testCanvas2 = new TCanvas("test2" ,"test2", 1600,800);
   testCanvas2->Divide(2,2);
   testCanvas2->cd(1);
   auto hTheta_cl = static_cast<TH1F*>(hTheta->Clone());
   hTheta_cl->Draw();
   hTheta_cl->SetTitle("a) Background Angle Distribution;angle [rad]; counts/1 mrad");
   hTheta_cl->SetTitleSize(hTheta_cl->GetTitleSize("x")*1.4, "xy");
   hTheta_cl->SetLabelSize(hTheta_cl->GetLabelSize("x")*1.25, "xy");

   auto pad2 = static_cast<TPad *>(testCanvas2->cd(2));
   pad2->SetLogy(1);
   auto hTheta2Log = static_cast<TH1F*>(hTheta2->Clone());
   hTheta2Log->SetTitle("b) LogY Photon Cherenkov Angle Distribution;angle [rad]; counts/1 mrad");
   hTheta2Log->SetTitleSize(hTheta_cl->GetTitleSize("x"), "xy");
   hTheta2Log->SetLabelSize(hTheta_cl->GetLabelSize("x"), "xy");
   hTheta2Log->Draw();


   auto pad3 = static_cast<TPad *>(testCanvas2->cd(3));
   //pad->SetLogy(1);
   auto hTheta2_cl = static_cast<TH1F*>(hTheta2->Clone());
   hTheta2_cl->SetTitleSize(hTheta_cl->GetTitleSize("x"), "xy");
   hTheta2_cl->SetLabelSize(hTheta_cl->GetLabelSize("x"), "xy");
   hTheta2_cl->Draw();
   hTheta2_cl->SetTitle("c) Photon Cherenkov Angle Distribution;angle [rad]; counts/1 mrad");
   auto pad4 = static_cast<TPad *>(testCanvas2->cd(4));
   //pad->SetLogy(1);
   auto hThetawg_cl = static_cast<TH1F*>(hThetawg->Clone());
   hThetawg_cl->SetTitleSize(hTheta_cl->GetTitleSize("x"), "xy");
   hThetawg_cl->SetLabelSize(hTheta_cl->GetLabelSize("x"), "xy");
   hThetawg_cl->SetTitle("d) (Weighted) Track Photon Cherenkov Angle Distribution;angle [rad]; counts/1 mrad");
   hThetawg_cl->Draw();
 }

 
 TFile *outFile = TFile::Open("background.root","RECREATE");
 
 outFile->WriteObject(hThetaRing,"hThetaRing");
 
 outFile->Write();
 outFile->Close();
 
}
//**********************************************************************************************************************************************************************************************************
Float_t GetFreonIndexOfRefraction(Float_t x)
  
{
  Float_t k = 1.177 + (0.0172)*x;
  return k;
}
//**********************************************************************************************************************************************************************************************************
Float_t GetQuartzIndexOfRefraction(Float_t x)
  
{
  Float_t k = TMath::Sqrt(1 + 46.411/(113.763556 - x) + 228.71/(328.51563 - x));
  return k;
}

//*********************************************************************************************************************************************************************************************************
Double_t BackgroundFunc(Double_t *x, Double_t *par)
{
 Double_t xx = x[0];
  
 Double_t f = par[0]*TMath::Tan(TMath::ASin(1.2903*TMath::Sin(xx)))*(1+TMath::Tan(TMath::ASin(1.2903*TMath::Sin(xx)))*TMath::Tan(TMath::ASin(1.2903*TMath::Sin(xx)))*1.2903*TMath::Cos(xx)/cos(asin(1.2903*TMath::Sin(xx))));
  
 return f;
}       
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


// get Ring-radius from Cherenkov-angle
double getRadiusFromCkov(double ckovAngle)
{

  //// refIndexFreon refIndexQuartz refIndexCH4
  double sin_ckov = static_cast<double>(TMath::Sin(ckovAngle));
  double sin_qz = static_cast<double>(sin_ckov*(refIndexFreon/refIndexQuartz));
  double sin_theta0 = static_cast<double>(sin_qz*(refIndexQuartz/refIndexCH4));

  double R_ckov = sin_ckov*(RadiatorWidth - EmissionLenght);
  double R_qz = sin_qz * QuartzWindowWidth;
  double R_0 = sin_theta0*CH4GapWidth;
  double R = static_cast<double>(R_ckov + R_qz + R_0);
  return R;
} 

