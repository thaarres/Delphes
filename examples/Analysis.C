//          ./Analysis delphes_output.root outfile.root



//remember to add your plots here
struct MyPlots
{
  TH1 *hJetBDeltaR;
  TH1 *hJetTrackDeltaR;
  TH1 *hJetEta;
  
  TH1 *hBJetEta;
  TH1 *hBJetPt;
  
  TH1 *hParticleEta;
  TH1 *hBQuarkEta;
  
  TH1 *hBJetIP_LeadingTrack;
  TH1 *hLightJetIP_LeadingTrack; 
  TH1 *hCJetIP_LeadingTrack;  
  
  TH1 *h_bJetTrackPt;
  TH1 *h_bJetLeadingTrackPt;
  TH1 *h_bJetTrackMultiplicity;
  
  TH1 *hBJetIP_first;
  TH1 *hBJetIP_second;
  TH1 *hBJetIP_third;
  
  TH1 *hCJetIP_first;
  TH1 *hCJetIP_second;
  TH1 *hCJetIP_third;
  
  TH1 *hLightJetIP_first;
  TH1 *hLightJetIP_second;
  TH1 *hLightJetIP_third;
  
  TH2 *h_IPversusTrackPt;
  
  
};

//------------------------------------------------------------------------------
class ExRootResult;
class ExRootTreeReader;
//------------------------------------------------------------------------------

//Function to compute dR
double DeltaR(TVector3 candidate1, TVector3 candidate2) 
{
  Double_t deta = candidate1.Eta()-candidate2.Eta();
  Double_t dphi = candidate1.Phi()-candidate2.Phi();
  return sqrt( deta*deta+dphi*dphi );
}
//------------------------------------------------------------------------------

//Function to compute the impact parameter
double ComputeIP(Track *track) {
  
  Double_t dxy;
  TVector3 vtrack;
  vtrack.SetPtEtaPhi(track->PT, track->Eta, track->Phi);

//    Double_t Bz = 3.8; //magnetic field
//    const Double_t c_light = 2.99792458E8; //c
//    
//    
//    Double_t r;
//    Double_t x_c, y_c, r_c, phi_c, phi_0, phi;
//    Double_t rcu, rc2, xd, yd;
//    Double_t X = track->X * 1.0E-3; //position in meters
//    Double_t Y = track->Y * 1.0E-3; //position in meters
// //   Double_t Z = track->Z * 1.0E-3; //position in meters
//    
//    r = vtrack.Pt() / (track->Charge * Bz) * 1.0E9/c_light; // helix radius in [m]
//    phi_0 = TMath::ATan2(vtrack.Py(), vtrack.Px()); // [rad] in [-pi, pi]
//   
//    // 2. helix axis coordinates
//    x_c = X + r*TMath::Sin(phi_0);
//    y_c = Y - r*TMath::Cos(phi_0);
//    r_c = TMath::Hypot(x_c, y_c);
//    phi_c = TMath::ATan2(y_c, x_c);
//    phi = phi_c;
//    if(x_c < 0.0) phi += TMath::Pi();
//    
//    rcu = TMath::Abs(r);
//    rc2 = r_c*r_c;
//    
//    // calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
//    xd = x_c*x_c*x_c - x_c*rcu*r_c + x_c*y_c*y_c;
//    xd = (rc2 > 0.0) ? xd / rc2 : -999;
//    yd = y_c*(-rcu*r_c + rc2);
//    yd = (rc2 > 0.0) ? yd / rc2 : -999;
//   //     zd = z + (TMath::Sqrt(xd*xd + yd*yd) - TMath::Sqrt(x*x + y*y))*pz/pt;
  
  // calculate impact paramater
  //dxy = (xd*vtrack.Py() - yd*vtrack.Px())/vtrack.Pt();
  //dxy = dxy * 1.0E2; //return value in cm
  
  dxy = (- track->X *vtrack.Py() + track->Y  *vtrack.Px())/track->PT;
  
  return dxy;
}

//------------------------------------------------------------------------------

//Function to compute the median

Double_t GetMedian(TH1 *h1) { 
   //compute the median for 1-d histogram h1 
   Int_t nbins = h1->GetXaxis()->GetNbins(); 
   Double_t *x = new Double_t[nbins]; 
   Double_t *y = new Double_t[nbins]; 
   for (Int_t i=0;i<nbins;i++) {
      x[i] = h1->GetXaxis()->GetBinCenter(i+1); 
      y[i] = h1->GetBinContent(i+1); 
   } 
   Double_t median = TMath::Median(nbins,x,y); 
   delete [] x; 
   delete [] y; 
   return median; 
}
//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, MyPlots *plots)
{
  
  THStack *stack;
  THStack *stack2;
  THStack *stack3;
  TLegend *legend;
  TLegend *legend2;
// TPaveText *comment;
  
  // AddHist1D(name, title, xlabel, ylabel, nxbins, xmin, xmax
  plots->hJetBDeltaR = result->AddHist1D(
    "hJetBDeltaR", "#DeltaR_{jet, b quark}",
    "#DeltaR_{jet, quark}", "number of jets",
      100, 0., 0.6);
  
  plots->hJetTrackDeltaR = result->AddHist1D(
    "hJetTrackDeltaR", "#DeltaR_{jet, track}",
    "R^{jet} - R^{track})", "number of jets",
     100, 0., 0.6);
  
  plots->hJetEta = result->AddHist1D(
    "jetEta", "Jet Eta",
    "Jet eta", "number of jets",
    100, -2.5, 2.5);
  
  plots->hBJetPt = result->AddHist1D(
    "hBJetPt", "b-jet p_{T}",
    "p_{T}", "number of jets",
    100, 20., 1000.);
    
      plots->hBJetEta = result->AddHist1D(
    "hBJetEta", "b-jet #eta",
    "#eta", "number of jets",
    100, -2.5, 2.5);
    
    plots->hBQuarkEta = result->AddHist1D(
    "hBQuarkEta", "b-quark #eta",
    "#eta", "number of b quarks",
    100, -2.5, 2.5);
  
  plots->hParticleEta = result->AddHist1D(
    "hParticleEta", "Gen Particle Eta",
    "Particle eta","number of jets",
    100, -2.5, 2.5);
  
  plots->hBJetIP_LeadingTrack = result->AddHist1D(
    "hBJetIP_LeadingTrack", "b-jet track IP",
    "b-jet track IP","number of tracks",
    400, -0.5, 0.5);
  
  plots->hCJetIP_LeadingTrack = result->AddHist1D(
    "hCJetIP_LeadingTrack", "c-jet track IP",
    "c-jet track IP","number of tracks",
    400, -0.5, 0.5);
  plots->hLightJetIP_LeadingTrack = result->AddHist1D(
    "hLightJetIP_LeadingTrack", "Light jet track IP",
    "Light jet track IP","number of tracks",
   400, -0.5, 0.5); 

plots->h_bJetTrackPt = result->AddHist1D(
    "h_bJetTrackPt", "b-jet track p_{T}",
    "Track p_{T}","number of tracks",
    60, 0., 60.);
    
    plots->h_bJetLeadingTrackPt = result->AddHist1D(
    "h_bJetLeadingTrackPt", "b-jet leading track p_{T}",
    "Leading track p_{T}","number of tracks",
    300, 0., 60.);   
    
    plots->h_bJetTrackMultiplicity = result->AddHist1D(
    "h_bJetTrackMultiplicity", "b-jet track multiplicity",
    "Track multiplicity","number of jets",
    300, 0., 25.); 
    
        plots->hLightJetIP_first = result->AddHist1D(
    "hLightJetIP_first", "Light-jet first track IP",
    "track IP","number of tracks",
       400, -0.5, 0.5);
      
          plots->hLightJetIP_second = result->AddHist1D(
    "hLightJetIP_second", "Light-jet second track IP",
    "track IP","number of tracks",
       400, -0.5, 0.5);
  
    plots->hLightJetIP_third = result->AddHist1D(
    "hLightJetIP_third", "Light-jet third track IP",
    "track IP","number of tracks",
     400, -0.5, 0.5);
      
          plots->hCJetIP_first = result->AddHist1D(
    "hCJetIP_first", "c-jet first track IP",
    "track IP","number of tracks",
        400, -0.5, 0.5);
      
          plots->hCJetIP_second = result->AddHist1D(
    "hCJetIP_second", "c-jet second track IP",
    "track IP","number of tracks",
        400, -0.5, 0.5);
  
    plots->hCJetIP_third = result->AddHist1D(
    "hCJetIP_third", "c-jet third track IP",
    "track IP","number of tracks",
      400, -0.5, 0.5);
      
    plots->hBJetIP_first = result->AddHist1D(
    "hBJetIP_first", "b-jet first track IP",
    "track IP","number of tracks",
      400, -0.5, 0.5);
    
    plots->hBJetIP_second = result->AddHist1D(
    "hBJetIP_second", "b-jet second track IP",
    "track IP","number of tracks",
       400, -0.5, 0.5);
  
    plots->hBJetIP_third = result->AddHist1D(
    "hBJetIP_third", "b-jet third track IP",
    "track IP","number of tracks",
     400, -0.5, 0.5);
     
     
     
    plots->h_IPversusTrackPt = result->AddHist2D(
    "h_IPversusTrackPt", "track IP versus p_{T}",
    "track IP","p_{T}",
     100, -0.3, 0.3, 100, 0.,100.);
     
  plots->hBJetIP_LeadingTrack->SetLineColor(kRed);
  plots->hCJetIP_LeadingTrack->SetLineColor(kBlue);
  plots->hLightJetIP_LeadingTrack->SetLineColor(kGreen);
  
  // book 1 stack of 2 histograms
  
  stack = result->AddHistStack("Leading track IP for b, c, light", "Leading track IP for b, c and light");
  stack->Add(plots->hBJetIP_LeadingTrack);
  stack->Add(plots->hLightJetIP_LeadingTrack);
  stack->Add(plots->hCJetIP_LeadingTrack);
  
  plots->hBJetIP_first->SetLineColor(kRed);
  plots->hBJetIP_second->SetLineColor(kBlue);
  plots->hBJetIP_third->SetLineColor(kGreen);
  
  plots->hCJetIP_first->SetLineColor(kBlue);
  plots->hLightJetIP_first->SetLineColor(kGreen);
  
  stack2 = result->AddHistStack("b-jet IP", "b-jet first, second and third highest IP track");
  stack2->Add(plots->hBJetIP_first);
  stack2->Add(plots->hBJetIP_second);
  stack2->Add(plots->hBJetIP_third);
  
  stack3 = result->AddHistStack("Highest IP", "Highest track IP for b, c and light ");
  stack3->Add(plots->hBJetIP_first);
  stack3->Add(plots->hCJetIP_first);
  stack3->Add(plots->hLightJetIP_first);
  
  legend2 = result->AddLegend(0.72, 0.86, 0.98, 0.98);
  legend2->AddEntry(plots->hBJetIP_first, "First", "f");
  legend2->AddEntry(plots->hBJetIP_second, "Second", "f");
  legend2->AddEntry(plots->hBJetIP_third, "Third", "f");
  legend2->Draw();
  result->Attach(stack2, legend2);
  
  // book legend for stack of 2 histograms
  
  legend = result->AddLegend(0.72, 0.86, 0.98, 0.98);
  legend->AddEntry(plots->hBJetIP_LeadingTrack, "b", "l");
  legend->AddEntry(plots->hCJetIP_LeadingTrack, "c ", "l");
  legend->AddEntry(plots->hLightJetIP_LeadingTrack, "light ", "l");
  
  
  // attach legend to stack (legend will be printed over stack in .eps file)
  
  result->Attach(stack, legend);
  
}

//------------------------------------------------------------------------------
void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
{
  
  
  
  //Branches
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  //   TClonesArray *branchEFlowTrack = treeReader->UseBranch("Track");
  
  //Objects 
  
  GenParticle *particle;
  
  Jet *jet;
  
  Track *track;
  Track *LeadingTrack;
  
  TVector3 vParticle; 
  TVector3 vJet;
  TVector3 vTrack;
  
  //Cuts
  Double_t pTMin  = 20;
  Double_t etaMin = 0.01;
  Double_t etaMax = 1.6;
  Double_t dRmin  = 0.5;
  Double_t PartondRmin  = 0.3;
  
  Double_t PartonPTMin = 1.0;
  Double_t PartonEtaMax  = 4.0;   
  
  
  Int_t i, j, k, pdgCode, pdgCodeMax, trackpT, trackpTMax; 
  
  Double_t LeadingTrackDxy; //Impact Parameter
  
  Long64_t allEntries = treeReader->GetEntries();
  Long64_t entry;  
  
  
  cout << "** Reading " << allEntries << " events" << endl;
  
  
  
  // Loop over events
  //     for(entry = 0; entry < 1; ++entry)
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    
    
    
    // Loop over all jets in event
    for(i = 0; i < branchJet->GetEntriesFast(); ++i)
    {
      Double_t PartonJetdR = -1;
      
      pdgCodeMax = -1;
      
      jet = (Jet*) branchJet->At(i);
      
      //Apply cuts
      if(jet->PT < pTMin || TMath::Abs(jet->Eta) > etaMax || TMath::Abs(jet->Eta) < etaMin ) continue;
      
      //Define a jet 3-vector
      vJet.SetPtEtaPhi(jet->PT, jet->Eta, jet->Phi);
      
      plots->hJetEta->Fill(jet->Eta);
      
      
      
      // Loop over GenParticles
      for(j = 0; j < branchParticle->GetEntriesFast(); ++j)
      {
	particle = (GenParticle*)branchParticle->At(j);
	
	if(particle->PT < PartonPTMin || TMath::Abs(particle->Eta) > PartonEtaMax) continue;
	
	vParticle.SetPtEtaPhi(particle->PT, particle->Eta, particle->Phi);
	
	plots->hParticleEta->Fill(particle->Eta);
	
	// Get jet flavor
	pdgCode = TMath::Abs(particle->PID);
	if(pdgCode == 5) plots->hBQuarkEta->Fill(particle->Eta);
	
	if(pdgCode != 21 && pdgCode > 5) continue;
	if(pdgCode == 21) pdgCode = 0;	  
	
	Double_t dR = DeltaR(vParticle, vJet);
	if(dR < PartondRmin)
	{
	  if(pdgCodeMax < pdgCode){
	    PartonJetdR = dR;
	    pdgCodeMax = pdgCode;
	  }
	}
      } 
      
      if(pdgCodeMax == 0) pdgCodeMax = 21;
      if(pdgCodeMax == -1) pdgCodeMax = 0;
      
      if( (pdgCodeMax != 21 && pdgCodeMax > 5) || pdgCodeMax == 0 ) continue;
      
      // If b jet
      if(pdgCodeMax == 5)
      {	  
	plots->hJetBDeltaR->Fill(PartonJetdR);
	plots->hBJetEta->Fill(jet->Eta);
	plots->hBJetPt->Fill(jet->PT);
	
      }
      
      trackpTMax = -1;
      Int_t NrTracks = 0;
      
      std::vector<double> IP;
      
      // Loop over tracks
      for(k = 0; k < branchEFlowTrack->GetEntriesFast(); ++k)
      {
	track = (Track*)branchEFlowTrack->At(k);
	if(track->PT < 1.0) continue;
	
	vTrack.SetPtEtaPhi(track->PT, track->Eta, track->Phi);
	Double_t dR = DeltaR(vTrack, vJet);
	
	// Get highest pT track in jet cone
	if(dR > dRmin)continue;
	
	Double_t dxy = ComputeIP(track);
	plots->h_IPversusTrackPt->Fill(dxy, track->PT);
	IP.push_back(TMath::Abs(dxy));
	
	NrTracks +=1;
	
	trackpT = track->PT;
	plots->hJetTrackDeltaR->Fill(dR);
	
	if(pdgCodeMax == 5){
		plots->h_bJetTrackPt->Fill(trackpT);
		}
		
	if(trackpTMax < trackpT){
	  trackpTMax = trackpT;
	  LeadingTrack = track;
	}
      }
      
      std::sort(IP.begin(),IP.end());
      std::reverse(IP.begin(),IP.end());
      
      
      // calculate Impact parameter
      LeadingTrackDxy = ComputeIP(LeadingTrack);
      
      // If b jet
      if(pdgCodeMax == 5)
      {	  
	plots->hBJetIP_LeadingTrack->Fill(LeadingTrackDxy);  
	plots->h_bJetLeadingTrackPt->Fill(LeadingTrack->PT);
	plots->h_bJetTrackMultiplicity->Fill(NrTracks);
	
	 if(IP.size()>0)plots->hBJetIP_first->Fill(IP[0]);  
      if(IP.size()>1)plots->hBJetIP_second->Fill(IP[1]);  
      if(IP.size()>2)plots->hBJetIP_third->Fill(IP[2]);  
      
      IP.clear();
      
      }
      
      // If c jet
      else if(pdgCodeMax == 4)
      {	  
	plots->hCJetIP_LeadingTrack->Fill(LeadingTrackDxy);
	
	if(IP.size()>0)plots->hCJetIP_first->Fill(IP[0]);  
      if(IP.size()>1)plots->hCJetIP_second->Fill(IP[1]);  
      if(IP.size()>2)plots->hCJetIP_third->Fill(IP[2]);  
      
      IP.clear();
      
      }
      
      // If light jet
      else if( (pdgCodeMax > 0) && (pdgCodeMax <=3 || pdgCodeMax ==21) )
      {	  
	plots->hLightJetIP_LeadingTrack->Fill(LeadingTrackDxy);
	
	if(IP.size()>0)plots->hLightJetIP_first->Fill(IP[0]);  
    if(IP.size()>1)plots->hLightJetIP_second->Fill(IP[1]);  
    if(IP.size()>2)plots->hLightJetIP_third->Fill(IP[2]);  
      
      IP.clear();
      
      }
      
    }	  
  }
  // plots->hBJetIP_LeadingTrack->Scale(1./plots->hBJetIP_LeadingTrack->Integral());
//   plots->hCJetIP_LeadingTrack->Scale(1./plots->hCJetIP_LeadingTrack->Integral());
//   plots->hLightJetIP_LeadingTrack->Scale(1./plots->hLightJetIP_LeadingTrack->Integral());
  
 //  plots->hLightJetIP_first->Scale(1./plots->hLightJetIP_first->Integral());
//   plots->hCJetIP_first->Scale(1./plots->hCJetIP_first->Integral());
//   plots->hBJetIP_first->Scale(1./plots->hBJetIP_first->Integral());
  
  Double_t median = GetMedian(plots->h_bJetTrackPt);
  std::cout<<"Median b-jetTrackPt"<<median<<std::endl;
  
  median = GetMedian(plots->hBJetPt);
  std::cout<<"Median hBJetPt"<<median<<std::endl;
  
  
  
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, MyPlots *plots)
{
  
  
  result->Print("png");
}

//------------------------------------------------------------------------------


void Analysis (const char *inputFile, const char *outputFile)
{
  gSystem->Load("libDelphes");
  
  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);
  chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_44_4to8k.root");
  chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_10_0to4k.root");
chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_10_12to16k.root");
chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_10_16to20k.root");
// chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_10_4to8k.root");
// chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_10_8to12k.root");
// chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_11_0to4k.root");
// chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_11_12to16k.root");
// chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_11_16to20k.root");
// chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_11_4to8k.root");
// chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_11_8to12k.root");
// chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_12_0to4k.root");
// chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_12_12to16k.root");
// chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_12_16to20k.root");
// chain->Add( "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/yangyong//data/Delphes/JetStudies_Phase_II_140PileUp_conf4/HH_bbWW/HH_bbWW_12_4to8k.root");

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();
  
  MyPlots *plots = new MyPlots;
  
  BookHistograms(result, plots);
  
  AnalyseEvents(treeReader, plots);
  
 // PrintHistograms(result, plots);
  
  result->Write(outputFile);
  
  cout << "** Exiting..." << endl;
  
  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}