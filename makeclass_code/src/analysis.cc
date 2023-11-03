#include "analysis.hh"


int main (int argc, char **argv)
{
  if (argc < 2)
    {
      printf("\nUsage: %s *.root\n\n",argv[0]);
      exit(0);
    }
  int const index = argc - 1;
  TString path[index];
  for ( int j = 0; j < argc-1; j++)
    {
      path[j] = (argv[j+1]);
      cout<<path[j]<<endl;
    }
  
  cout<<"File name = "<<path[0]<<endl;
  
  TString file_details= path[0];
  std::string str_input(file_details.Data());
  //str_input += "/*.root";
  std::string str(file_details.Data());
  std::replace(str.begin(), str.end(), '_', ' ');  // replace '_' by ' '

  vector<string> array;
  stringstream ss(str);
  string temp;
  while (ss >> temp)
    array.push_back(temp); 

  TChain *chain = new TChain("Events");
  /*
   //Data
   for (int k = 1; k < 474 ; k++){     //0.49/fb of era B
       std::string str_input_file(str_input + "/output_"+ std::to_string(k)+"_Friend.root");
       if(gSystem->AccessPathName(str_input_file.c_str())==0){
           chain -> Add(str_input_file.c_str()); // Add the root file to the TChain chain
           cout << " adding " << str_input_file << endl;
           }
       }
  */ 
  for (int k = 0; k < argc-1 ; k++)
    {
  
  chain -> Add(path[k]); // Add the root file to the TChain chain
  cout << " adding " << path[k] << endl;
   
    }

  p = new Events(chain);

  int numev = p->fChain->GetEntries();
  cout << "Total number of events to be processed = " << numev << endl;

  get_event(0) ;

  for( int i = 0 ; i < numev ; i ++ )
     {
     get_event(i);
     cout << "Event " << i << endl; 
     /* 
     for(int genpart=0;genpart<p->nGenPart;genpart++){
         //if(p->GenPart_status[genpart]==1){
             int pdgid_original = p->GenPart_pdgId[genpart];
             pdgid_original_distribution->Fill(pdgid_original); 
 
             int index = p->GenPart_genPartIdxMother[genpart];
             //cout << "Mother index = " << index << endl;
             if(abs(p->GenPart_pdgId[index])==4900113 && abs(p->GenPart_eta[index])<2.4){ //omega
             //if((abs(p->GenPart_pdgId[index])==999999) && abs(p->GenPart_eta[index])<2.4){
                 ROOT::Math::XYZVector genpart_vertex;
                 genpart_vertex.SetXYZ(p->GenPart_vertex_x[genpart],p->GenPart_vertex_y[genpart],p->GenPart_vertex_z[genpart]);
                 ROOT::Math::XYZVector mother_vertex;
                 mother_vertex.SetXYZ(p->GenPart_vertex_x[index],p->GenPart_vertex_y[index],p->GenPart_vertex_z[index]);
                 ROOT::Math::XYZVector displacement = genpart_vertex - mother_vertex;

                 int pdgid = p->GenPart_pdgId[genpart];
                 pdgid_distribution->Fill(pdgid);   

                 Double_t d3d = displacement.R();
                 ROOT::Math::PtEtaPhiMVector mother_P(p->GenPart_pt[index], p->GenPart_eta[index], p->GenPart_phi[index], p->GenPart_mass[index]);
                 Double_t E = mother_P.E();
                 Double_t P = mother_P.P();
                 Double_t gamma = E/p->GenPart_mass[index];
                 Double_t beta = P/E;
                 Double_t ctau = d3d/(beta*gamma);
                
                 cout << "Event " << i << endl;        
  
                 cout << "Parent mass = " << p->GenPart_mass[index] << endl;
                 
                 cout << "Vertex d3d = " << d3d << endl;
                 cout << "E = " << E << endl;
                 cout << "gamma = " << gamma << endl;
                 cout << "beta = " << beta << endl;
                 cout << "ctau = " << ctau << endl;
               
                 ctau_distribution->Fill(ctau);
                 d3d_distribution->Fill(d3d);

                 }   
             //}
         }
     //gen part loop ends
     */
    
     for(int genpart=0;genpart<p->nGenPart;genpart++){
         if(abs(p->GenPart_pdgId[genpart])==13 && abs(p->GenPart_eta[genpart])<2.4){ //gen muons
             gen_muon_pT->Fill(p->GenPart_pt[genpart]);
             gen_muon_vertex_dxy->Fill(sqrt(p->GenPart_vertex_x[genpart]*p->GenPart_vertex_x[genpart] + p->GenPart_vertex_y[genpart]*p->GenPart_vertex_y[genpart]));
         } 
  
     }

     /*
     //match reco muon to gen muon
     for(int recomu=0; recomu<p->nMuon; recomu++){
       int genpartidx = p->Muon_genPartIdx[recomu]; 
       int genpart = p->nGenPart;
       if(genpartidx > genpart - 1) continue;
       if(abs(p->GenPart_pdgId[genpartidx])==13){
          cout << "Reco muon dxy = " << abs(p->Muon_dxy[recomu]) << endl; 
       }
     }
     //gen part loop ends
     */
     
     vector<TVector3> genpart_vertices0;
     vector<TVector3> genpart_vertices;
     vector<TVector3> genpart_vertices1;

     genpart_vertices0.clear();
     genpart_vertices.clear();
     genpart_vertices1.clear();
 
     for(int genpart0=0;genpart0<(int)p->nGenPart;genpart0++){  
      TVector3 genpart_vertex0; 
      genpart_vertex0.SetXYZ(p->GenPart_vertex_x[genpart0],p->GenPart_vertex_y[genpart0],p->GenPart_vertex_z[genpart0]);
      genpart_vertices0.push_back(genpart_vertex0);
     }
  
     for(int genpart=0;genpart<(int)p->nGenPart;genpart++){
	    int leadingmupt_index=-1, lowestmupt_index=-1;
      int recomupt_index1=-1, recomupt_index2=-1;
      int leadingrecomupt_index=-1, lowestrecomupt_index=-1;
      int recomuplus_index = -1, recomuminus_index = -1;
      int muplus_index = -1, muminus_index = -1;    
        
      if(p->GenPart_status[genpart]!=1 || fabs(p->GenPart_pdgId[genpart])!=13 || fabs(p->GenPart_pt[genpart])<=3.0 || fabs(p->GenPart_eta[genpart])>=2.4) continue;
 
      int mu1_match = 0; int mu2_match = 0;

      TVector3 genpart_vertex;
	    genpart_vertex.SetXYZ(p->GenPart_vertex_x[genpart],p->GenPart_vertex_y[genpart],p->GenPart_vertex_z[genpart]);

      for(int recomu=0;recomu<(int)p->nMuon;recomu++){
	     if(p->Muon_genPartIdx[recomu]==genpart){ 
          mu1_match++;
          recomupt_index1 = recomu;
          }
	    }
      int match=0;
	    for(int genpart1=genpart+1;genpart1<(int)p->nGenPart;genpart1++){
	     if(p->GenPart_status[genpart1]!=1 || fabs(p->GenPart_pdgId[genpart1])!=13 || fabs(p->GenPart_pt[genpart1])<=3.0 || fabs(p->GenPart_eta[genpart1])>=2.4) continue;
       TVector3 genpart_vertex1;
	     for(int recomu=0;recomu<(int)p->nMuon;recomu++){
	       if(p->Muon_genPartIdx[recomu]==genpart1){
              mu2_match++;
              recomupt_index2 = recomu;
         }
	      }
        genpart_vertex1.SetXYZ(p->GenPart_vertex_x[genpart1],p->GenPart_vertex_y[genpart1],p->GenPart_vertex_z[genpart1]);
          //int genpartidxmother = p->GenPart_genPartIdxMother[genpart];
          //int genpartidxmother1 = p->GenPart_genPartIdxMother[genpart1]; 
	      if(genpart_vertex==genpart_vertex1){
            //int genpart_mother_match = 0;
            //for(int genpart2=0;genpart2<(int)p->nGenPart;genpart2++){
              //if(p->GenPart_genPartIdxMother[genpart2]==p->GenPart_genPartIdxMother[genpart1]){genpart_mother_match++;}
            //}
            int check_vertex = 0;
            for(int index=0;index<genpart_vertices0.size(); index++){
              if(genpart_vertices0.at(index)==genpart_vertex1){check_vertex++;}
            }

            if(check_vertex==2){
	            match++;
              leadingmupt_index=genpart; lowestmupt_index=genpart1;
              leadingrecomupt_index = recomupt_index1; lowestrecomupt_index = recomupt_index2;
              
              if(p->GenPart_pt[genpart1]>p->GenPart_pt[genpart]){ leadingmupt_index=genpart1; lowestmupt_index=genpart; leadingrecomupt_index = recomupt_index2; lowestrecomupt_index = recomupt_index1;}

              recomuplus_index = recomupt_index1; recomuminus_index = recomupt_index2; muplus_index = genpart; muminus_index = genpart1;
              if(p->Muon_charge[recomupt_index2]== 1){recomuplus_index = recomupt_index2; recomuminus_index = recomupt_index1; muplus_index = genpart1; muminus_index = genpart;}
            } 
	       }

       }
      
       if(match>0){ 

	      int check_previous_vertex=0;
	      for(int index=0;index<genpart_vertices.size(); index++){
	        if(genpart_vertices.at(index)==genpart_vertex){
              check_previous_vertex++;
          }
        }
        if(check_previous_vertex!=0){cout << "check_previous_vertex= " << check_previous_vertex << endl;} 
	      if(check_previous_vertex==0){
	       genpart_vertices.push_back(genpart_vertex);
         if(mu1_match>0 && mu2_match>0 && leadingmupt_index != -1 && lowestmupt_index != -1 && leadingrecomupt_index != -1 && lowestrecomupt_index != -1){
	       //cout<<"Matched both reco muons for this event"<<endl;
	       genpart_vertices1.push_back(genpart_vertex);
	       }
        }
       }   
     } //Gen particle loop ends

     //Check matching between generated SV and min chi2 SV
     for(int genpart=0;genpart<(int)genpart_vertices1.size();genpart++){
	     int match_dR=0, match_dR_distance=0, match_dR_distance1=0;
       int match_dR_distance_darkqcd=0;
       TVector3 genpart_vertex = genpart_vertices1.at(genpart);
 
       for(int sv=0; sv<(int)p->nSV; sv++){
         TVector3 sv_vertex;
	       sv_vertex.SetXYZ(p->SV_x[sv],p->SV_y[sv],p->SV_z[sv]);
         double dR_distance = (genpart_vertex-sv_vertex).Perp();
         if(dR_distance < 0.1){
           gen_dxy_reco_dxy_SV->Fill(genpart_vertex.Perp(), sv_vertex.Perp());
           gen_d3d_reco_d3d_SV->Fill(genpart_vertex.Mag(), sv_vertex.Mag());
         }

       }
     }

       
     for(int muonsv=0; muonsv < p->nmuonSV; muonsv++){
        int mu1_idx = p->muonSV_mu1index[muonsv];
        int mu2_idx = p->muonSV_mu2index[muonsv];

        int nmuon = p->nMuon;
        cout << "Event " << i << endl;
        cout << "mu1_idx = " << mu1_idx << endl;
        cout << "mu2_idx = " << mu2_idx << endl;
        cout << "nMuon = " << nmuon << endl;

        if(abs(mu1_idx) > nmuon - 1 || abs(mu2_idx) > nmuon - 1) continue;
 
        int mu1_genpartidx = p->Muon_genPartIdx[mu1_idx];
        int mu2_genpartidx = p->Muon_genPartIdx[mu2_idx];

        int ngenpart = p->nGenPart;
        cout << "mu1_genpartidx = " << mu1_genpartidx << endl;
        cout << "mu2_genpartidx = " << mu2_genpartidx << endl;
        cout << "ngenpart = " << ngenpart << endl;

        if(abs(mu1_genpartidx) > ngenpart - 1 || abs(mu2_genpartidx) > ngenpart - 1) continue; 

        ROOT::Math::XYZVector genpart_muon_vertex;
        genpart_muon_vertex.SetXYZ(p->GenPart_vertex_x[mu1_genpartidx],p->GenPart_vertex_y[mu1_genpartidx],p->GenPart_vertex_z[mu1_genpartidx]);
        ROOT::Math::XYZVector genpart_muon_vertex1;
        genpart_muon_vertex1.SetXYZ(p->GenPart_vertex_x[mu2_genpartidx],p->GenPart_vertex_y[mu2_genpartidx],p->GenPart_vertex_z[mu2_genpartidx]);  

        ROOT::Math::XYZVector reco_muon_vertex;
        reco_muon_vertex.SetXYZ(p->muonSV_x[muonsv], p->muonSV_y[muonsv], p->muonSV_z[muonsv]);        
        /*
        ROOT::Math::XYZVector reco_muon_vertex_orig;
        reco_muon_vertex_orig.SetXYZ(p->muonSV_x[muonsv], p->muonSV_y[muonsv], p->muonSV_z[muonsv]);

        ROOT::Math::XYZVector reco_muon_vertex;
        int match = 0;

        for(int sv=0; sv < p->nSV; sv++){
           ROOT::Math::XYZVector SV;
           SV.SetXYZ(p->SV_x[sv], p->SV_y[sv], p->SV_z[sv]);
           if((SV - reco_muon_vertex_orig).Rho() < 0.1){ 
              reco_muon_vertex.SetXYZ(p->SV_x[sv], p->SV_y[sv], p->SV_z[sv]);
              match ++;
           }
        }

        cout << "match = " << match;

        if(match == 0) continue;
        */
 
        int mu1_genpartidx_mother = p->GenPart_genPartIdxMother[mu1_genpartidx];
        int mu2_genpartidx_mother = p->GenPart_genPartIdxMother[mu2_genpartidx];

        if(mu1_genpartidx_mother > ngenpart - 1 || mu2_genpartidx_mother > ngenpart - 1) continue;

        if(genpart_muon_vertex == genpart_muon_vertex1 && abs(p->GenPart_pdgId[mu1_genpartidx])==13 && abs(p->GenPart_pdgId[mu2_genpartidx])==13){
        //if(genpart_muon_vertex == genpart_muon_vertex1 && mu1_genpartidx_mother == mu2_genpartidx_mother && p->GenPart_pdgId[mu1_genpartidx_mother] == 999999){

           cout << "Gen vertex d3d = " << genpart_muon_vertex.R() << endl;
           cout << "Reco muonSV d3d = " << reco_muon_vertex.R() << endl;
           cout << "Gen vertex dxy = " << genpart_muon_vertex.Rho() << endl;
           cout << "Reco muonSV dxy = " << reco_muon_vertex.Rho() << endl;

           genmuonsv_d3d_distribution->Fill(genpart_muon_vertex.R());
           recomuonsv_d3d_distribution->Fill(reco_muon_vertex.R());
           genminusreco_d3d_distribution->Fill(genpart_muon_vertex.R() - reco_muon_vertex.R());
           gen_d3d_reco_d3d->Fill(genpart_muon_vertex.R(), reco_muon_vertex.R());

           genmuonsv_dxy_distribution->Fill(genpart_muon_vertex.Rho());
           recomuonsv_dxy_distribution->Fill(reco_muon_vertex.Rho());
           genminusreco_dxy_distribution->Fill(genpart_muon_vertex.Rho() - reco_muon_vertex.Rho());
           gen_dxy_reco_dxy->Fill(genpart_muon_vertex.Rho(), reco_muon_vertex.Rho());

           if((genpart_muon_vertex.Rho() - reco_muon_vertex.Rho()) > 4.0){
               cout << "nMuon = " << p->nMuon << endl;
           }
           }
        

        } 
     //muonsv loop ends
     
     }

    /*
    TFile f("output/output_scenarioA_mpi_10_mA_3p33_ctau_10_confirm.root","recreate"); 
    //TFile f("output/output_scenarioA_mpi_1_mA_0p33_ctau_0p1_confirm.root","recreate"); 
    //TFile f("output/output_vector_2_10_1_1_confirm.root","recreate");
    ctau_distribution->Write();
    d3d_distribution->Write();
    pdgid_original_distribution->Write();
    pdgid_distribution->Write();
    

    f.Close();
    */
      
    TFile f("output/ul_pu_new/outputSV_UL_scenarioA_mpi_4_mA_1p33_ctau_10.root","recreate");
    //TFile f("output/outputSV_scenarioA_mpi_4_mA_1p33_ctau_10_with_confirm.root","recreate");
    //ctau_distribution->Write();
    //d3d_distribution->Write();
    //pdgid_original_distribution->Write();
    //pdgid_distribution->Write();
    genmuonsv_d3d_distribution->Write();
    recomuonsv_d3d_distribution->Write();
    genmuonsv_dxy_distribution->Write();
    recomuonsv_dxy_distribution->Write();
    genminusreco_d3d_distribution->Write();
    genminusreco_dxy_distribution->Write();
    gen_d3d_reco_d3d->Write();
    gen_dxy_reco_dxy->Write();
    gen_muon_pT->Write();
    gen_muon_vertex_dxy->Write();
    gen_d3d_reco_d3d_SV->Write();
    gen_dxy_reco_dxy_SV->Write();
    
  
    return 0;
 }
    








