//-------------------------------------------Functions------------------------------------------------------------------//    

void formatCorrHist(TH2D* hist) {

    hist->GetXaxis()->SetTitle("#Delta#eta");
    hist->GetYaxis()->SetTitle("#Delta#phi");
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleOffset(1.3); 
   
    return;
    
}

void formatCorrHist(TH1D* hist) {

    hist->GetXaxis()->SetTitle("#Delta#phi");
    hist->GetXaxis()->SetTitleOffset(1.3);
   
   return;
    
}


void formatCorrHist(TH2D* hist, TString title) {

    hist->GetXaxis()->SetTitle("#Delta#eta");
    hist->GetYaxis()->SetTitle("#Delta#phi");
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleOffset(1.3); 
    hist->GetXaxis()->SetTitle("#Delta#eta");
    hist->GetYaxis()->SetTitle("#Delta#phi");
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleOffset(1.3); 
    hist->SetNameTitle(title, title);
    
    return;
    
}








    /*for(int j = 0; j < 2; j++){ //begin block to obtain the US and LS histograms and calculate the Norm factor from side bands
    
        oneDHistos[j] = (TH1D*)file->Get(oneDStrings[j]);
        
        TAxis *xAxis = oneDHistos[j]->GetXaxis();
        Int_t binMassLow = xAxis->FindBin(sideBandLow);            //get normalization factors to scale the LS distribution
        Int_t binMassHigh = xAxis->FindBin(sideBandHigh);
        
        integralSideBandCounts[j] = oneDHistos[j]->Integral(binMassLow, binMassHigh);
        
        Int_t binMassLow = xAxis->FindBin(massLow);            
        Int_t binMassHigh = xAxis->FindBin(massHigh);
        
        integralCounts[j] = oneDHistos[j]->Integral(binMassLow, binMassHigh);
        
        cout << oneDStrings[j] << "   SideBandCounts: " << integralSideBandCounts[j] << endl;  
                
    
        oneDHistos[j]->GetXaxis()->SetTitle("Invariant Mass GeV/c^{2}");
        oneDHistos[j]->GetYaxis()->SetTitle("counts");
        oneDHistos[j]->GetYaxis()->SetTitleOffset(1.2);
        oneDHistos[j]->Sumw2();
        oneDHistos[j]->SetMarkerStyle(20);
        oneDHistos[j]->Fit(g1, "R0");
        oneDHistos[j]->Fit(e1, "R0+");
        g1->GetParameters(&par[0]);
        e1->GetParameters(&par[3]);
        
        fun->SetParameters(par);
        oneDHistos[j]->Fit(fun, "R+");
        //oneDHistos[1]->SetMarkerColor(2);
        oneDHistos[j]->Draw();
        //oneDHistos[1]->Draw("same");
        
        str1 = path + oneDStrings[j] + fileType;
        cout << str1 << endl;
        c->SaveAs(str1);
    
    }  //end block to obtain the US and LS histograms and calculate the Norm factor from side bands  
    
    
    
    cout << "Low and high bins in mass range: " << binMassLow << "   " << binMassHigh << endl;
    cout << oneDStrings[0] << "   SignalCounts: " << integralCounts[0] << endl; 
    cout << oneDStrings[1] << "   LSBGCounts: " << integralCounts[1] << endl; 
    cout << "Significance (with unscaled LS): " << (integralCounts[0]-integralCounts[1])/(TMath::Sqrt(integralCounts[0])) << endl;
    
    
    oneDHistos[2]->Add(oneDHistos[0], oneDHistos[1], 1, -1);     //US - unscaled LS
    oneDHistos[2]->Draw();
    str1 = path + oneDStrings[2] + fileType;                     
    cout << str1 << endl;
    c->SaveAs(str1);
    
    
    oneDHistos[3] = (TH1D*)file->Get(oneDStrings[1]);            //This creates the histogram to store the scaled LS
    
    LSScaleFactor = integralSideBandCounts[0]/integralSideBandCounts[1];
    
    cout << "LS Scale Factor: " << LSScaleFactor << endl;
    
    oneDHistos[3]->Scale(integralSideBandCounts[0]/integralSideBandCounts[1]);   //scale the LS histogram
    oneDHistos[3]->SetTitle(oneDStrings[3]);
    oneDHistos[3]->GetXaxis()->SetTitle("Invariant Mass GeV/c^{2}");
    oneDHistos[3]->GetYaxis()->SetTitle("counts");
    oneDHistos[3]->GetYaxis()->SetTitleOffset(1.2);
    
    
    
    oneDHistos[3]->Draw();
    str1 = path + oneDStrings[3] + fileType;                     //store the LS histogram as a .png
    cout << str1 << endl;
    c->SaveAs(str1);
     
    //Deepa stuff-------------------------------
    
    oneDHistos[1]->SetMarkerColor(2);
    oneDHistos[0]->Draw();
    oneDHistos[3]->Draw("same");
    
    //oneDHistos[0]->Rebin();
   // oneDHistos[3]->Rebin();
    
    str1 = path + "Deepa" + fileType;
    cout << str1 << endl;
    c->SaveAs(str1);
    
    //Deepa stuff---------------------------------
    
    TAxis *xAxis = oneDHistos[0]->GetXaxis();
    Int_t binMassLow = xAxis->FindBin(massLow);            
    Int_t binMassHigh = xAxis->FindBin(massHigh);
    
    integralCounts[0] = oneDHistos[0]->Integral(binMassLow, binMassHigh); 
    
    TAxis *xAxis = oneDHistos[3]->GetXaxis();
    Int_t binMassLow = xAxis->FindBin(massLow);            
    Int_t binMassHigh = xAxis->FindBin(massHigh);
    
    integralCounts[1] = oneDHistos[3]->Integral(binMassLow, binMassHigh); 
    
    oneDHistos[4]->Add(oneDHistos[0], oneDHistos[3], 1, -1);
    oneDHistos[4]->SetTitle(oneDStrings[4]);
    oneDHistos[4]->GetXaxis()->SetTitle("Invariant Mass GeV/c^{2}");
    oneDHistos[4]->GetYaxis()->SetTitle("counts");
    oneDHistos[4]->GetYaxis()->SetTitleOffset(1.2);
    
    //FITTING////////////
    oneDHistos[4]->Fit(g1, "R0");
    oneDHistos[4]->Fit(e1, "R0+");
    g1->GetParameters(&par[0]);
    e1->GetParameters(&par[3]);
        
    fun->SetParameters(par);
    oneDHistos[4]->Fit(fun, "R+");
    ////////////
    
    oneDHistos[4]->Draw();
    str1 = path + oneDStrings[4] + fileType;
    
    cout << str1 << endl;
    
    c->SaveAs(str1);
    
    cout << oneDStrings[0] << "   SignalCounts: " << integralCounts[0] << endl; 
    cout << oneDStrings[3] << "   ScaledLSBGCounts: " << integralCounts[1] << endl; 
    
    cout << "Significance (with scaled LS): " << (integralCounts[0]-integralCounts[1])/(TMath::Sqrt(integralCounts[0])) << endl;    
    
    //file->Close();
    
    for(int i = 0; i < 5; i++){
    
        oneDHistos[i]->Write();
    }

//------------------------------------------
//END INVARIANT MASS INFORMATION    
//------------------------------------------*/