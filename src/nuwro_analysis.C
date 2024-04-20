// Selection of CCQE in NuWro
// The momentum of hit neutron is less than 0.3 GeV
// Reconstruction of neutrino energy from the neutron kinematics

struct EveInfo{
  double nMom, pMom , TrueEnu, RecoEnu ;
};

EveInfo kinematic(event *ev)
{
    EveInfo info;

    double E_true = ev->in[0].E(); //True Energy of neutrino

    double El = 0.0; 
    double lx = 0.0; 
    double ly = 0.0; 
    double lz = 0.0;         // lepton E and P
    double Ep = 0.0;
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;         // Proton E and P

    for (int n = 0 ; n < ev->post.size(); n++)    // particle loop
    {
        if (ev->post[n].proton())   // for proton
        {
            Ep = ev->post[n].E();
            px = ev->post[n].p().x;
            py = ev->post[n].p().y;
            pz = ev->post[n].p().z;
        }

        if (ev->post[n].lepton())  // for muon
        {
            El = ev->post[n].E();
            lx = ev->post[n].p().x;
            ly = ev->post[n].p().y;
            lz = ev->post[n].p().z;
        }
    }

    // calculating MX and  MY
    // for argon
    // int Np = 18; int Nn = 22;                                //number of proton and neutron 
    // double mp = 938.27208816; double mn = 939.5654205;       //mass of proton and neutron
    // double BE = 343.84; double e = 24.78;                    // binding energy and separation energy
    // double MX = (Np * mp) + (Nn * mn) - BE;                  // calculating M_A
    // double MY = MX - mn + e;                                 // M_(A-1)

    // // for carbon
    int Np = 6; int Nn = 6;                              //number of proton and neutron 
    double mp = 938.27208816; double mn = 939.5654205;   //mass of proton and neutron
    double BE = 92.16; double e = 27.139;                // binding energy and separation energy
    double MX = (Np * mp) + (Nn * mn) - BE;              // calculating M_A
    double MY = MX - mn + e;  

    // mag of proton momentum
    double pmom = sqrt(pow(px,2) + pow(py,2) + pow(pz,2));

    // recostruction of neutron momentum
    double Sq_nT = pow((lx+px),2) + pow((ly+py),2); // (perpendicular comp. of n (in paper pT))
    double A     = MX + lz + pz - El - Ep ;
    double nL    = (pow(A,2) - Sq_nT - pow(MY,2))/(2*A);

    // total momentum
    double Pn = std::sqrt(Sq_nT + pow(nL,2));

    // neutrino energy reconstruction
    double Enu = lz + pz - nL;



    info.nMom    = Pn;
    info.pMom    = pmom;
    info.TrueEnu = E_true;
    info.RecoEnu = Enu;

    return info;
}

void analysis()
{
    TFile *f = new TFile("ND_C.nuwro.root");  //input file
    TTree *t = (TTree*)f->Get("treeout");

    TFile *output = new TFile("ND_C.nuwro.out.root","RECREATE"); //output file
    
    TTree *tree1 = new TTree("neutron", "neutron");
    TTree *tree2 = new TTree("selection","selection");

    double pnFull ; double pnReco; int pID; int label; double Ereco; double Etrue; int sID; int EveID; double ppFull;

    tree1->Branch("Process",&pID);
    tree1->Branch("MomN",&pnFull);
    tree1->Branch("MomP",&ppFull);

    tree2->Branch("EventID",&EveID);
    tree2->Branch("MomN",&pnReco);
    tree2->Branch("label",&label);
    tree2->Branch("ScatteringID",&sID);
    tree2->Branch("RecoEnu",&Ereco);
    tree2->Branch("TrueEnu",&Etrue);

    event *e = new event();
    t->SetBranchAddress("e",&e);

    for (int i = 0; i < t->GetEntries(); i++)
    {
        t->GetEntry(i);

        int no_proton = 0; int no_pion = 0;   // no of p and pi

        for (int k =0 ; k < e->post.size(); k++)
        {
            if (e->post[k].proton()){no_proton++;}
            if (e->post[k].pion()){no_pion++;}
        }

        if (no_pion ==0 && no_proton ==1)
        {
            auto result = kinematic(e);

            pID    = e->dyn;
            pnFull = (result.nMom)/1000;
            ppFull = (result.pMom)/1000;

            if ((result.nMom < 300) && (result.pMom > 200))
            {
                // signal

                if (e->flag.cc && e->flag.qel)
                {
                    EveID  =i;
                    label  = 1;
                    pnReco = (result.nMom)/1000;
                    sID    = e->dyn;
                    Ereco  = (result.RecoEnu)/1000;
                    Etrue  = (result.TrueEnu)/1000;
                }

                // background
                else
                {
                    EveID  =i;
                    label  = 0;
                    pnReco = (result.nMom)/1000;
                    sID    = e->dyn;
                    Ereco  = (result.RecoEnu)/1000;
                    Etrue  = (result.TrueEnu)/1000; 
                }

                tree2->Fill();
            }

            tree1->Fill();
        }
    }     // end of event loop

    tree1->Write();
    tree2->Write();

    output->Write();
    output->Close();

    cout <<"*****************DONE***************"<<endl;

}