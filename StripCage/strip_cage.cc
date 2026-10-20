#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TLine.h>
#include <TColor.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include "TSystem.h"
#include "Garfield/ComponentElmer.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/Random.hh"
#include "Garfield/MediumGas.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/DriftLineRKF.hh"
#include <csignal>
#include <csetjmp>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TEllipse.h>
#include <TBox.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TPad.h>


// OpenMP intentionally not used: parallel AvalancheMC drifts caused
// cross-thread corruption. All modes run serial.




using namespace Garfield;






void PlotFieldPointMap(Garfield::Sensor& sensor,
                       const std::string& out_dir,
                       const std::string& tag,
                       double d    = 0.1,      // characteristic spacing [cm]
                       double xlo  = -2.0, double xhi =  2.0,
                       double ylo  = -2.0, double yhi =  2.0,
                       double zmin = 0.0,  double zmax = 0.0,   // 0,0 = autoscale
                       double zprobe = 0.005,
                       int canvasW = 700) {


  const double r = 0.375 * d;                 // radius -> diameter = (3/4)*d


  std::vector<double> px, py, pe;
  double emin =  1e30, emax = -1e30;








  for (double x = xlo + 0.5 * d; x <= xhi; x += d) {
    for (double y = ylo + 0.5 * d; y <= yhi; y += d) {
      double ex = 0., ey = 0., ez = 0., v = 0.;
      Garfield::Medium* med = nullptr;
      int status = 0;
      sensor.ElectricField(x, y, zprobe, ex, ey, ez, v, med, status);
      if (status != 0 || med == nullptr) continue;   // inside metal or off-mesh
      const double e = std::sqrt(ex * ex + ey * ey + ez * ez);
      px.push_back(x); py.push_back(y); pe.push_back(e);
      if (e < emin) emin = e;
      if (e > emax) emax = e;
    }
  }
  if (px.empty()) {
    std::cout << "PlotFieldPointMap: no valid points sampled.\n";
    return;
  }
  if (zmax > zmin) { emin = zmin; emax = zmax; }        // manual clip
  if (emax <= emin) emax = emin + 1.0;




  std::cout << "  field point map: " << px.size() << " points, |E| range "
            << emin << " .. " << emax << " V/cm\n";






  gStyle->SetPalette(kBird);
  const int nCol = TColor::GetNumberOfColors();


  // keep the dots round
  // canvas aspect follows the data aspect
  int canvasH = int(canvasW * (yhi - ylo) / (xhi - xlo));
  if (canvasH < 300)  canvasH = 300;
  if (canvasH > 2400) canvasH = 2400;


  TCanvas* c = new TCanvas(("cFieldPts_" + tag).c_str(),
                           "E field point map", canvasW, canvasH);
  c->SetRightMargin(0.18);
  c->SetLeftMargin(0.14);


  // supplies the axes and the palette axis on the right.
  TH2D* hF = new TH2D(("hFieldFrame_" + tag).c_str(),
                      ";x [cm];y [cm]", 10, xlo, xhi, 10, ylo, yhi);
  hF->SetStats(0);
  hF->SetBinContent(1, 1, emin);      // one entry so ROOT draws the palette
  hF->SetMinimum(emin);
  hF->SetMaximum(emax);
  hF->GetZaxis()->SetTitle("|E| [V/cm]");
  hF->Draw("COLZ");
  gPad->Update();


  // wipe the dummy bin, keep frame + palette
  TBox* blank = new TBox(xlo, ylo, xhi, yhi);
  blank->SetFillColor(kWhite);
  blank->Draw();


  for (size_t i = 0; i < px.size(); ++i) {
    int ic = int(nCol * (pe[i] - emin) / (emax - emin));
    if (ic < 0)      ic = 0;
    if (ic >= nCol)  ic = nCol - 1;
    const int col = TColor::GetColorPalette(ic);
    TEllipse* el = new TEllipse(px[i], py[i], r, r);
    el->SetFillColor(col);
    el->SetLineColor(col);
    el->Draw();
  }




  // strip planes for reference
  TLine* sL = new TLine(-1.5, ylo, -1.5, yhi);
  TLine* sR = new TLine( 1.5, ylo,  1.5, yhi);
  sL->SetLineColor(kGray + 2); sL->SetLineWidth(2); sL->Draw();
  sR->SetLineColor(kGray + 2); sR->SetLineWidth(2); sR->Draw();




  gPad->RedrawAxis();
  c->SaveAs((out_dir + "/strip_field_point_map_" + tag + ".pdf").c_str());
  delete c;
}












int main(int argc, char * argv[]) {
    TApplication app("app", &argc, argv);
    gROOT->SetBatch(true);
















    // Create timestamped output directory
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time), "%Y-%m-%d_%H-%M-%S");
    std::string out_dir = "Outputs/" + ss.str();
   
    gSystem->MakeDirectory("Outputs");
    gSystem->MakeDirectory(out_dir.c_str());
    std::cout << "\nAll outputs will be saved to: " << out_dir << "\n";
    std::cout << "Running single-threaded (OpenMP disabled)\n";


    MediumGas gas;


    ComponentElmer elm("strip_cage/mesh.header", "strip_cage/mesh.elements", "strip_cage/mesh.nodes",
                       "dielectrics.dat", "strip_cage/strip_cage.result", "cm");
    elm.SetMedium(0, &gas);
    elm.EnablePeriodicityZ();




    Sensor sensor;
    sensor.AddComponent(&elm);
   
    // The strips now have rounded corners, eliminating the mathematical
    // singularities. We can now safely set the sensor area all the way out to x = +/- 2.0,
    // exactly matching field_cage.cc, allowing electrons to drift safely behind the strips!
    sensor.SetArea(-1.5, -16.0, -100.0, 1.5, 16.0, 100);


    AvalancheMC drift;
    drift.SetSensor(&sensor);








    drift.SetDistanceSteps(0.02);








    drift.SetTimeWindow(0.0, 500000.0);
      
    drift.EnableAttachment();     // ← ADD (was missing)
    drift.EnableDiffusion();      // ← ADD
    ViewDrift driftView;
    driftView.SetArea(-2.0, -16.0, -1.0, 2.0, 16.0, 1.0);






    // ======================================================================
    // SWITCH: 0 = Visual Diagnostic,  1 = Multi-Gas Spatial Sweep,
    //         2 = Diffusion test,     3 = Sweep + readout response,
    //         4 = Ripple test (D_T vs launch x0),
    //         6 = Full E-field dump,  7 = Micro E-field probe (3x3 grid)
    // ======================================================================
    int switch_mode = 1;


    if (switch_mode == 0) {
        std::cout << "\n=========================================\n";
        std::cout << "MODE 0: VISUAL DIAGNOSTIC (67% CF4)\n";
        std::cout << "=========================================\n";








        const std::string gas_file_m0 = "gas_tables/cf4_80_iso_20.gas";
        if (!gas.LoadGasFile(gas_file_m0)) {
            std::cerr << "FATAL: could not load " << gas_file_m0 << "\n";
            gSystem->Exit(1);
        }
        std::cout << "-> Gas: " << gas_file_m0 << "\n";
        drift.EnablePlotting(&driftView);








        const int nElectrons = 500;
        int leaked = 0;
       
        std::cout << "Drifting " << nElectrons << " random electrons for visualization...\n";
        for (int e = 0; e < nElectrons; ++e) {
            double x0 = -1.4 + RndmUniform()*2.8;
            double y0 = -14 + RndmUniform()*28;
            double z0 = 0.005, t0 = 0.0;




            std::cout << "  [e=" << e << "] x0=" << x0 << " y0=" << y0 << std::flush;




   
                        drift.DriftElectron(x0, y0, z0, t0);
                        const int nPoints = drift.GetNumberOfDriftLinePoints();
                        if (nPoints > 0) {
                            double x1, y1, z1, t1;
                            drift.GetDriftLinePoint(nPoints - 1, x1, y1, z1, t1);
                            if ((y1 < -15.98) && (std::abs(x1) < 1.5)) leaked++;
                            std::cout << " end (" << x1 << ", " << y1 << ") t=" << t1 << " ns, npts=" << nPoints << "\n";




                        }
                        std::cout << " -> ok\n";
                   






           
        }
        std::cout << "RESULTS: " << leaked << " out of " << nElectrons << " leaked past the field cage!\n\n";






        // 1. Drift Lines
        TCanvas* c1 = new TCanvas("c1", "Electron Drift Lines", 600, 600);
        driftView.SetCanvas(c1);
        driftView.Plot(true);
        TLine* stripL = new TLine(-1.5, -16.0, -1.5, 16.0);
        TLine* stripR = new TLine( 1.5, -16.0,  1.5, 16.0);
        stripL->SetLineColor(kGray+2); stripL->SetLineWidth(2); stripL->Draw();
        stripR->SetLineColor(kGray+2); stripR->SetLineWidth(2); stripR->Draw();
        c1->SaveAs((out_dir + "/strip_drift_lines_67_33.pdf").c_str());
        std::cout << "-> Saved " << out_dir << "/strip_drift_lines_67_33.pdf\n";








        // 2. Field Map
        TCanvas* cMap = new TCanvas("cMap", "Electric Field Contours", 800, 800);
        ViewField fieldView;
        fieldView.SetSensor(&sensor);
        fieldView.SetCanvas(cMap);
        fieldView.SetPlane(0, 0, 1, 0, 0, 0.005);
        fieldView.SetArea(-1.6, -16, 1.6, -14);
        fieldView.SetNumberOfContours(50);
        fieldView.SetNumberOfSamples2d(400, 400);
        //fieldView.SetVoltageRange(-1002., -980.);
        fieldView.PlotContour("v");
        TLine* stripL2 = new TLine(-1.5, -16.0, -1.5, 16.0);
        TLine* stripR2 = new TLine( 1.5, -16.0,  1.5, 16.0);
        stripL2->SetLineColor(kGray+2); stripL2->SetLineWidth(2); stripL2->Draw();
        stripR2->SetLineColor(kGray+2); stripR2->SetLineWidth(2); stripR2->Draw();
        cMap->SaveAs((out_dir + "/strip_field_map_contours_67_33.pdf").c_str());
        std::cout << "-> Saved " << out_dir << "/strip_field_map_contours_67_33.pdf\n";








            // Zoomed: drift region near the middle, clipped colour range so the
        // bulk uniformity is actually visible.
        PlotFieldPointMap(sensor, out_dir, "zoom", 0.05,
                        -1.9, 1.9, -2.0, 2.0, 25.0, 45.0);








        // Top of the cage -- this is the region driving the efficiency upturn.
        PlotFieldPointMap(sensor, out_dir, "topend", 0.05,
                        -1.9, 1.9, 12.0, 16.0, 25.0, 45.0);








        // Full height, coarse grid, tall canvas.
        PlotFieldPointMap(sensor, out_dir, "full", 0.25,
                        -1.9, 1.9, -15.9, 15.9, 25.0, 45.0, 0.005, 320);










        // 3. Drift Velocity Graph
        ViewMedium mediumView;
        mediumView.SetMedium(&gas);
        TCanvas* cVel = new TCanvas("cVel", "Drift Velocity", 800, 600);
        mediumView.SetCanvas(cVel);
        mediumView.PlotElectronVelocity('e');
        cVel->SaveAs((out_dir + "/strip_drift_velocity_67_33.pdf").c_str());
        std::cout << "-> Saved " << out_dir << "/strip_drift_velocity_67_33.pdf\n";
















        // 4. Mirror Symmetry Diagnostic
        // 4. Mirror Symmetry Diagnostic — scan the full drift volume
        double sum_diff = 0.0, sum_mag = 0.0;
        int nPts = 0;
        for (double yy = -15.9; yy <= 15.9; yy += 0.2) {
            for (double xx = 0.05; xx <= 1.49; xx += 0.05) {
                double ex1, ey1, ez1, v1, ex2, ey2, ez2, v2;
                Medium *m1 = nullptr, *m2 = nullptr;
                int s1 = 0, s2 = 0;
                elm.ElectricField( xx, yy, 0.005, ex1, ey1, ez1, v1, m1, s1);
                elm.ElectricField(-xx, yy, 0.005, ex2, ey2, ez2, v2, m2, s2);
                if (s1 != 0 || s2 != 0) continue;      // skip points outside the mesh
                sum_diff += std::abs(v1 - v2);
                sum_mag  += (std::abs(v1) + std::abs(v2)) / 2.0;
                ++nPts;
            }
        }
        std::cout << "-> Mirror Symmetry Diagnostic: "
                  << 100.0 * (1.0 - (sum_diff / sum_mag))
                  << "%  (" << nPts << " point pairs)\n";
















        // ==========================================================
        // 5. Electric Field Uniformity Diagnostic
        // ==========================================================
        std::cout << "\n--- FIELD UNIFORMITY DIAGNOSTIC ---\n";
       
        double sum_e = 0.0;
        double sum_e2 = 0.0;
        int n_uniform_pts = 0;
       
        // Scan parameters
        double x_step = 0.1; // Check a column every 1 mm
        double y_step = 0.1; // Check down the Y axis every 2 mm
        double z_val = 0.005;






        for (double x = -1.3; x <= 1.301; x += x_step) {
            for (double y = -15.0; y <= 15.0; y += y_step) {
               
                double ex = 0., ey = 0., ez = 0., v = 0.;
                Medium* m = nullptr;
                int status = 0;
               
                elm.ElectricField(x, y, z_val, ex, ey, ez, v, m, status);
               
                if (status == 0) { // status == 0 means the point is safely inside the gas mesh
                    double emag = std::sqrt(ex*ex + ey*ey + ez*ez);
                    sum_e += emag;
                    sum_e2 += (emag * emag);
                    n_uniform_pts++;
                }
            }
        }




        if (n_uniform_pts > 1) {
            double avg_e = sum_e / n_uniform_pts;
           
            // Standard deviation formula: sqrt( (Sum(x^2) - (Sum(x)^2 / N)) / (N - 1) )
            double variance = (sum_e2 - ((sum_e * sum_e) / n_uniform_pts)) / (n_uniform_pts - 1);
            if (variance < 0.0) variance = 0.0; // Prevent NaN from floating point rounding errors
            double std_dev = std::sqrt(variance);
















            std::cout << "Sampled " << n_uniform_pts << " valid points inside X:[-1.49, 1.49], Y:[-15.0, 15.0]\n";
            std::cout << "Average Field (|E|) : " << avg_e << " V/cm\n";
            std::cout << "Standard Deviation  : " << std_dev << " V/cm\n";
           
            // Optional: Print it as a percentage of the average field (Coefficient of Variation)
            std::cout << "Field Fluctuation   : " << (std_dev / avg_e) * 100.0 << " %\n";
        } else {
            std::cout << "Error: Not enough valid mesh points found in the scanned region.\n";
        }
        std::cout << "-----------------------------------\n\n";
















     
        {
        double ex = 0., ey = 0., ez = 0., v = 0.;
        Medium* m = nullptr;
        int status = 0;
        elm.ElectricField(1.0, -15.0, 0.005, ex, ey, ez, v, m, status);
        std::cout << "EDGE PROBE: |E|=" << sqrt(ex*ex+ey*ey+ez*ez) << " V/cm status=" << status << "\n";
        }
































    } else if (switch_mode == 1) {
        std::cout << "\n=========================================\n";
        std::cout << "MODE 1:  MULTI-GAS SWEEP\n";
        std::cout << "=========================================\n";




        // ==========================================================
        // 1. CALCULATE AVERAGE ELECTRIC FIELD
        // ==========================================================
        double sum_e = 0.0;
        int n_uniform_pts = 0;
        for (double x = -1.490; x <= 1.491; x += 0.1) {
            for (double y = -15.0; y <= 15.0; y += 0.2) {
                double ex = 0., ey = 0., ez = 0., v = 0.;
                Medium* m = nullptr;
                int status = 0;
                elm.ElectricField(x, y, 0.005, ex, ey, ez, v, m, status);
               
                if (status == 0) { // status == 0 means valid gas region
                    sum_e += std::sqrt(ex*ex + ey*ey + ez*ez);
                    n_uniform_pts++;
                }
            }
        }
        // Fallback to 33.0 if no points found, otherwise use calculated average
        double avg_e = (n_uniform_pts > 0) ? (sum_e / n_uniform_pts) : 33.0;
        std::cout << "-> Calculated Average Field for Sweep: " << avg_e << " V/cm\n";
































        // ==========================================================
        // 2. CONFIGURE GRAPHS
        // ==========================================================
        std::vector<std::string> o2_tags = {
            "cf4_67_iso_33_o2_0p000pct", "cf4_67_iso_33_o2_0p010pct",
            "cf4_67_iso_33_o2_0p100pct", "cf4_67_iso_33_o2_1p000pct",
            "cf4_67_iso_33_o2_5p000pct"};
        std::vector<std::string> o2_labels = {
            "0ppm", "100ppm", "1000ppm", "1pct", "5pct"};
        double y_step = 0.1;
        // Sweep extent. Change these and the loop below follows.
        const double Y_MIN = -15.0;
        const double Y_MAX =  15.0;
        int colors[] = {kBlack, kBlue, kGreen+2, kOrange+1, kRed};
















        // MultiGraphs for Spatial Efficiency
        TMultiGraph* mg = new TMultiGraph();
        mg->SetTitle("Total Spatial Efficiency vs Spawn Y;Spawn Y Position [cm];Electrons Passed [%]");
       
        TMultiGraph* mgCurves = new TMultiGraph();
        mgCurves->SetTitle("Spatial Efficiency Fits;Spawn Y Position [cm];Electrons Passed [%]");
       
        auto legend = new TLegend(0.80, 0.52, 0.985, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.035);
















        // Transverse Diffusion Graph (Using dynamic average field)
        TGraph* grDiff = new TGraph();
        grDiff->SetTitle(Form("Transverse Diffusion at %.1f V/cm;CF4 Concentration [%%];Transverse Diffusion [cm^{1/2}]", avg_e));
        grDiff->SetMarkerStyle(20); grDiff->SetMarkerColor(kRed); grDiff->SetLineColor(kRed); grDiff->SetLineWidth(2);






        double eta = 0.0;
        gas.ElectronAttachment(0.0, avg_e, 0.0, 0., 0., 0., eta);
        std::cout << "attachment eta = " << eta << " /cm";
        if (eta > 1e-9) std::cout << "   (lambda = " << 1.0 / eta << " cm)";
        std::cout << "\n";












        // NEW: Drift Velocity Graph (Using dynamic average field)
        TGraph* grVel = new TGraph();
        grVel->SetTitle(Form("Drift velocity vs cf4/isobutane concentration at %.1f V/cm;CF4 Concentration [%%];Drift Velocity [cm/ns]", avg_e));
        grVel->SetMarkerStyle(20); grVel->SetMarkerColor(kMagenta+2); grVel->SetLineColor(kMagenta+2); grVel->SetLineWidth(2);
















        int graph_idx = 0;
















        // Single drift tracker for the serial sweep.
        AvalancheMC local_drift;
        local_drift.SetSensor(&sensor);
        local_drift.SetDistanceSteps(0.02);
        local_drift.SetTimeWindow(0.0, 500000.0);   // ns — adjust to your real drift time
        local_drift.EnableAttachment();
        local_drift.EnableDiffusion();
















        // ==========================================================
        // 3. RUN GAS SWEEP
        // ==========================================================
        for (int i = 0; i < (int)o2_tags.size(); ++i) {
            const std::string tag   = o2_tags[i];
            const std::string label = o2_labels[i];
            std::string gas_file = "gas_tables_v2/" + tag + ".gas";
            if (!gas.LoadGasFile(gas_file)) {
                std::cerr << "FATAL: could not load " << gas_file << "\n";
                gSystem->Exit(1);
            }
















            std::string out_name = out_dir + "/epass_" + tag + ".txt";            
            std::ofstream outFile(out_name);
            outFile << "Spawn_Y\tElectrons_Passed_Pct\n";
            std::string lost_name = out_dir + "/lost_" + tag + ".txt";            
            std::ofstream lostFile(lost_name);
            lostFile << "Spawn_Y\tx0\tx_end\ty_end\tt_end\tnpts\n";
            std::cout << "\nTesting Gas: 67/33 CF4/iC4H10 + O2 " << label << "\n";















            // Record Transverse Diffusion at avg_e
            double dl = 0, dt = 0;
            gas.ElectronDiffusion(0.0, avg_e, 0.0, 0., 0., 0., dl, dt);
            grDiff->SetPoint(graph_idx, i, dt);
















            // NEW: Record Drift Velocity at avg_e
            double vx = 0., vy = 0., vz = 0.;
            gas.ElectronVelocity(0.0, avg_e, 0.0, 0., 0., 0., vx, vy, vz);
            double v_mag = std::sqrt(vx*vx + vy*vy + vz*vz);
            grVel->SetPoint(graph_idx, i, v_mag);
















            TGraph* grSpatial = new TGraph();
            int pt_idx = 0;
































           
            const int nY = int(std::lround((Y_MAX - Y_MIN) / y_step)) + 1;
            for (int iy = 0; iy < nY; ++iy) {
                const double y_spawn = Y_MIN + y_step * iy;
















                const int nElectrons = 500;
                int leaked = 0;
















                for (int e = 0; e < nElectrons; ++e) {
                    const double x0 = -1.4 + (2.8 * e) / (nElectrons - 1.0);
                    local_drift.DriftElectron(x0, y_spawn, 0.005, 0.0);

                    const int nPoints = local_drift.GetNumberOfDriftLinePoints();
                    double x1 = 0, y1 = 0, z1 = 0, t1 = 0;
                    if (nPoints > 0) local_drift.GetDriftLinePoint(nPoints - 1, x1, y1, z1, t1);

                    if (nPoints > 0 && (y1 < -15.98) && (std::abs(x1) < 1.5)) {
                        leaked++;
                    } else {
                        lostFile << y_spawn << "\t" << x0 << "\t" << x1 << "\t"
                                 << y1 << "\t" << t1 << "\t" << nPoints << "\n";
                    }
                }
























                double pass_pct = 100.0 * leaked / nElectrons;
                outFile << y_spawn << "\t" << pass_pct << "\n";
                grSpatial->SetPoint(pt_idx++, y_spawn, pass_pct);
                std::cout << "  Y = " << y_spawn << " -> " << pass_pct
                          << "% passed\n" << std::flush;
            }
















            outFile.close();
















            grSpatial->SetLineColor(colors[graph_idx]);
            grSpatial->SetLineWidth(2);
            grSpatial->SetMarkerStyle(20);
            grSpatial->SetMarkerSize(0.25);      // quarter of the ROOT default
            grSpatial->SetMarkerColor(colors[graph_idx]);
            mg->Add(grSpatial, "LP");
            legend->AddEntry(grSpatial, Form("O2 %s", label.c_str()), "lp");
















            // INDIVIDUAL FITS & PLOTS
            TGraph* grForFit = (TGraph*)grSpatial->Clone(Form("grFit_%d", i));
            TF1* fitFunc = new TF1(Form("fit_%d", i), "[0]/(x - [1]) + [2]", -15.0, 15.0);
            fitFunc->SetParameters(10.0, -20.0, 50.0);
            fitFunc->SetLineColor(colors[graph_idx]);
            grForFit->Fit(fitFunc, "Q");
















            TCanvas* cSingle = new TCanvas(Form("cSingle_%d", i), Form("O2 %s Efficiency", label.c_str()), 800, 600);
            cSingle->SetGrid();
            grForFit->SetTitle(Form("Spatial Efficiency: O2 %s;Spawn Y Position [cm];Electrons Passed [%%]", label.c_str()));
            grForFit->Draw("AP");
            cSingle->SaveAs((out_dir + "/" + tag + "_EPass_v_Y.pdf").c_str());
            delete cSingle;
















            TGraph* grCurveOnly = new TGraph(fitFunc);
            grCurveOnly->SetLineColor(colors[graph_idx]);
            grCurveOnly->SetLineWidth(2);
            mgCurves->Add(grCurveOnly, "L");
           
            graph_idx++;
        }
















        // ==========================================================
        // 4. FINAL OUTPUT PLOTS
        // ==========================================================
















        // Plot Spatial Efficiency TOTAL
        TCanvas* cEff = new TCanvas("cEff", "Total Spatial Efficiency", 1000, 600);
        cEff->SetRightMargin(0.22);
        cEff->SetGrid();
        mg->SetMinimum(0.0);
        mg->SetMaximum(100.0);
        mg->Draw("A");
       
        mg->GetXaxis()->SetLimits(Y_MIN, Y_MAX);
        mg->GetXaxis()->SetNdivisions(510);
        gPad->Modified();
        gPad->Update();
        legend->Draw();
        cEff->SaveAs((out_dir + "/Total_EPass_v_Y.pdf").c_str());
















        // Plot Spatial Efficiency TOTAL CURVES
        TCanvas* cCurves = new TCanvas("cCurves", "Total Fit Curves", 1000, 600);
        cCurves->SetRightMargin(0.22);
        cCurves->SetGrid();
        mgCurves->SetMinimum(0.0);
        mgCurves->SetMaximum(100.0);
        mgCurves->Draw("A");
        mgCurves->GetXaxis()->SetLimits(Y_MIN, Y_MAX);
        mgCurves->GetXaxis()->SetNdivisions(510);
        gPad->Modified();
        gPad->Update();
        legend->Draw();
        cCurves->SaveAs((out_dir + "/TotalCurves_EPass_v_Y.pdf").c_str());
















        // Plot Transverse Diffusion
        TCanvas* cDiff = new TCanvas("cDiff", "Transverse Diffusion", 800, 600);
        cDiff->SetGrid();
        grDiff->Draw("APL");
        cDiff->SaveAs((out_dir + "/strip_transverse_diffusion_vs_cf4.pdf").c_str());
















        // NEW: Plot Drift Velocity
        TCanvas* cVelSweep = new TCanvas("cVelSweep", "Drift Velocity vs CF4", 800, 600);
        cVelSweep->SetGrid();
        grVel->Draw("APL");
        cVelSweep->SaveAs((out_dir + "/strip_drift_velocity_vs_cf4.pdf").c_str());
        std::cout << "-> Saved " << out_dir << "/strip_drift_velocity_vs_cf4.pdf\n";
    } else if (switch_mode == 2) {
        // ==================================================================
        //  MODE 2 - TRANSVERSE / LONGITUDINAL DIFFUSION MEASUREMENT
        //
        //  Launch N electrons from a single point on the axis and follow the
        //  cloud down the cage, sampling the beam width at fixed heights.
        //
        //  The key idea: do NOT just look at where the electrons land. Each
        //  drift line already contains the full trajectory, so one set of N
        //  electrons gives sigma(L) at EVERY depth for free. Endpoint-only
        //  would need a separate run per drift length.
        //
        //  Diffusion theory says   sigma_x = D_T * sqrt(L)
        //  so sigma_x plotted against sqrt(L) is a straight line through the
        //  origin whose SLOPE is D_T in cm^(1/2). That slope, compared with
        //  what MediumGas tabulates, validates the whole transport chain.
        // ==================================================================
        // Loop over the SAME five gases the mode-1 sweep uses, so the
        // measured diffusion is reported for every O2 level we tested.
        std::vector<std::string> o2_tags = {
            "cf4_67_iso_33_o2_0p000pct", "cf4_67_iso_33_o2_0p010pct",
            "cf4_67_iso_33_o2_0p100pct", "cf4_67_iso_33_o2_1p000pct",
            "cf4_67_iso_33_o2_5p000pct"};
        std::vector<std::string> o2_labels = {"0ppm","100ppm","1000ppm","1pct","5pct"};
        std::vector<double> o2_pct = {0.0, 0.01, 0.10, 1.0, 5.0};   // O2 % for the x-axis




        const int    nElectrons = 5000;    // 10x the old 500 for a clean D_T fit
                                           // (statistics use all; only nPlot are drawn)
        const double x_start    =  0.0;    // on axis
        const double y_start    =  0.0;    // exact middle of the cage
        const double z_start    =  0.005;
        const double y_stop     = -15.0;   // deepest sampling plane
        const double dy_check   =  0.5;    // sampling plane spacing [cm]




        // Average |E| over the drift volume. The field is set by Elmer and is
        // gas-independent, so compute it once, outside the gas loop.
        double sum_e2 = 0.0;
        int n_pts2 = 0;
        for (double x = -1.490; x <= 1.491; x += 0.1) {
            for (double y = -15.0; y <= 15.0; y += 0.2) {
                double ex = 0., ey = 0., ez = 0., v = 0.;
                Medium* m = nullptr;
                int status = 0;
                elm.ElectricField(x, y, 0.005, ex, ey, ez, v, m, status);
                if (status == 0) {
                    sum_e2 += std::sqrt(ex * ex + ey * ey + ez * ez);
                    ++n_pts2;
                }
            }
        }
        const double avg_e = (n_pts2 > 0) ? (sum_e2 / n_pts2) : 33.0;
        std::cout << "-> Average field for mode 2: " << avg_e << " V/cm\n";




        const int nCheck = int(std::lround((y_start - y_stop) / dy_check));
        std::vector<double> yCheck(nCheck);
        for (int c = 0; c < nCheck; ++c) yCheck[c] = y_start - dy_check * (c + 1);

        const int nPlot = 100;   // draw only the first nPlot drift lines (first gas only)




        // Summary table across all gases: measured vs tabulated D_T and D_L.
        std::ofstream fsum(out_dir + "/diffusion_vs_o2.txt");
        fsum << "O2_label\tO2_pct\tE_Vcm\tDT_meas\tDT_meas_err\tDT_tab\tDT_ratio"
                "\tDL_meas\tDL_tab\tDL_ratio\tv_drift_cmns\n";

        // Graphs: measured and tabulated D_T (and D_L) vs O2 %.
        TGraph* grDTmeas = new TGraph();  TGraph* grDTtab = new TGraph();
        TGraph* grDLmeas = new TGraph();  TGraph* grDLtab = new TGraph();

        for (int ig = 0; ig < (int)o2_tags.size(); ++ig) {
            const std::string gas_file = "gas_tables_v2/" + o2_tags[ig] + ".gas";
            if (!gas.LoadGasFile(gas_file)) {
                std::cerr << "FATAL: could not load " << gas_file << "\n";
                gSystem->Exit(1);
            }
            std::cout << "\n==================================================\n";
            std::cout << "GAS " << (ig+1) << "/" << o2_tags.size()
                      << " : 67/33 CF4/iC4H10 + O2 " << o2_labels[ig]
                      << "  (" << gas_file << ")\n";

            // Per-gas accumulators (reset each gas).
            std::vector<long long> nAt(nCheck, 0);
            std::vector<double> sX(nCheck, 0.0), sX2(nCheck, 0.0);
            std::vector<double> sT(nCheck, 0.0), sT2(nCheck, 0.0);
            std::vector<double> sX4(nCheck, 0.0);

            AvalancheMC dif;
            dif.SetSensor(&sensor);
            dif.SetDistanceSteps(0.01);
            dif.SetTimeWindow(0.0, 500000.0);
            dif.EnableDiffusion();
            // Attachment deliberately OFF: this mode measures the pure diffusion
            // coefficient (cloud spread). Enabling attachment would preferentially
            // remove longer-path electrons, biasing the surviving-sample sigma_x.
            // Adding O2 barely changes D_T itself; its loss is captured in mode 1.
            if (ig == 0) dif.EnablePlotting(&driftView);   // only draw lines for the first gas

            TH1D* hEnd = new TH1D(Form("hEnd_%d", ig),
                                  "Final x after full drift;x [cm];Electrons", 150, -1.5, 1.5);

            int nComplete = 0, nAbsorbed = 0;
            for (int e = 0; e < nElectrons; ++e) {
                if (ig == 0 && e == nPlot) dif.DisablePlotting();
                dif.DriftElectron(x_start, y_start, z_start, 0.0);
                const int nPts = dif.GetNumberOfDriftLinePoints();
                if (nPts < 2) continue;

                double xE, yE, zE, tE;
                dif.GetDriftLinePoint(nPts - 1, xE, yE, zE, tE);
                if (yE < y_stop) { ++nComplete; hEnd->Fill(xE); }
                else             { ++nAbsorbed; }

                // Walk the trajectory once, crossing sampling planes in order.
                int c = 0;
                double xP, yP, zP, tP;
                dif.GetDriftLinePoint(0, xP, yP, zP, tP);
                for (int p = 1; p < nPts && c < nCheck; ++p) {
                    double x, y, z, t;
                    dif.GetDriftLinePoint(p, x, y, z, t);
                    while (c < nCheck && yP > yCheck[c] && y <= yCheck[c]) {
                        const double dyseg = yP - y;
                        const double f = (dyseg > 1e-12) ? (yP - yCheck[c]) / dyseg : 0.0;
                        const double xi = xP + f * (x - xP);
                        const double ti = tP + f * (t - tP);
                        nAt[c] += 1;
                        sX[c]  += xi;   sX2[c] += xi * xi;   sX4[c] += xi * xi * xi * xi;
                        sT[c]  += ti;   sT2[c] += ti * ti;
                        ++c;
                    }
                    xP = x; yP = y; zP = z; tP = t;
                }
                if ((e + 1) % 100 == 0)
                    std::cout << "  " << (e + 1) << " / " << nElectrons << "\r" << std::flush;
            }
            std::cout << "\n  Reached the bottom: " << nComplete
                      << "   absorbed early: " << nAbsorbed << "\n";

            // ---- tabulated (pure-gas) prediction for THIS gas ----
            double dl_tab = 0.0, dt_tab = 0.0;
            gas.ElectronDiffusion(0.0, avg_e, 0.0, 0., 0., 0., dl_tab, dt_tab);
            double vx_t = 0., vy_t = 0., vz_t = 0.;
            gas.ElectronVelocity(0.0, avg_e, 0.0, 0., 0., 0., vx_t, vy_t, vz_t);
            const double v_tab = std::sqrt(vx_t * vx_t + vy_t * vy_t + vz_t * vz_t);

            // Per-gas depth profile.
            std::ofstream df(out_dir + "/diffusion_profile_" + o2_tags[ig] + ".txt");
            df << "L_cm\tn\tmean_x\tsigma_x\tDT_local\tmean_t_ns\tsigma_t_ns\tDL_local\texcess_kurt\n";

            TGraph* grSig  = new TGraph();
            TGraph* grDTL  = new TGraph();
            TGraph* grSigT = new TGraph();
            int k = 0;
            for (int c = 0; c < nCheck; ++c) {
                if (nAt[c] < 50) continue;
                const double N  = double(nAt[c]);
                const double L  = y_start - yCheck[c];
                const double mx = sX[c] / N;
                double vx = sX2[c] / N - mx * mx;            if (vx < 0) vx = 0;
                const double sx = std::sqrt(vx);
                const double mt = sT[c] / N;
                double vt = sT2[c] / N - mt * mt;            if (vt < 0) vt = 0;
                const double st = std::sqrt(vt);
                const double kurt = (vx > 0) ? (sX4[c] / N) / (vx * vx) - 3.0 : 0.0;
                const double sL   = st * v_tab;
                const double DT_l = sx / std::sqrt(L);
                const double DL_l = sL / std::sqrt(L);
                df << L << "\t" << nAt[c] << "\t" << mx << "\t" << sx << "\t" << DT_l
                   << "\t" << mt << "\t" << st << "\t" << DL_l << "\t" << kurt << "\n";
                grSig ->SetPoint(k, std::sqrt(L), sx);
                grSigT->SetPoint(k, std::sqrt(L), sL);
                grDTL ->SetPoint(k, L, DT_l);
                ++k;
            }
            df.close();

            // ---- fits: sigma = D * sqrt(L) through the origin ----
            TF1* fT = new TF1(Form("fT_%d",ig), "[0]*x", 0.0, std::sqrt(y_start - y_stop) * 1.05);
            fT->SetParameter(0, dt_tab); grSig->Fit(fT, "Q");
            const double DT_fit = fT->GetParameter(0);
            const double DT_err = fT->GetParError(0);

            TF1* fL = new TF1(Form("fL_%d",ig), "[0]*x", 0.0, std::sqrt(y_start - y_stop) * 1.05);
            fL->SetParameter(0, dl_tab); grSigT->Fit(fL, "Q");
            const double DL_fit = fL->GetParameter(0);

            const double DT_ratio = (dt_tab > 0 ? DT_fit / dt_tab : 0);
            const double DL_ratio = (dl_tab > 0 ? DL_fit / dl_tab : 0);

            std::cout << "  ---- DIFFUSION (E = " << avg_e << " V/cm, O2 " << o2_labels[ig] << ") ----\n";
            std::cout << "    measured D_T = " << DT_fit << " +/- " << DT_err
                      << " | tabulated D_T = " << dt_tab << " | ratio " << DT_ratio << "\n";
            std::cout << "    measured D_L = " << DL_fit
                      << " | tabulated D_L = " << dl_tab << " | ratio " << DL_ratio << "\n";
            std::cout << "    drift velocity (table) = " << v_tab << " cm/ns\n";

            fsum << o2_labels[ig] << "\t" << o2_pct[ig] << "\t" << avg_e << "\t"
                 << DT_fit << "\t" << DT_err << "\t" << dt_tab << "\t" << DT_ratio << "\t"
                 << DL_fit << "\t" << dl_tab << "\t" << DL_ratio << "\t" << v_tab << "\n";

            grDTmeas->SetPoint(ig, o2_pct[ig], DT_fit);
            grDTtab ->SetPoint(ig, o2_pct[ig], dt_tab);
            grDLmeas->SetPoint(ig, o2_pct[ig], DL_fit);
            grDLtab ->SetPoint(ig, o2_pct[ig], dl_tab);

            // Per-gas sigma-vs-sqrt(L) plot.
            TCanvas* cS = new TCanvas(Form("cS_%d",ig), "sigma vs sqrt(L)", 800, 600);
            cS->SetGrid();
            grSig->SetTitle(Form("Transverse spread, O2 %s;#sqrt{L} [cm^{1/2}];#sigma_{x} [cm]", o2_labels[ig].c_str()));
            grSig->SetMarkerStyle(20); grSig->SetMarkerSize(0.5);
            grSig->SetMarkerColor(kRed + 1); grSig->SetLineColor(kRed + 1);
            grSig->Draw("AP");
            cS->SaveAs((out_dir + "/diffusion_sigma_vs_sqrtL_" + o2_tags[ig] + ".pdf").c_str());
            delete cS;

            // Drift lines only for the first gas.
            if (ig == 0) {
                TCanvas* cDr = new TCanvas("cDr", "Diffusion Drift Lines", 600, 600);
                driftView.SetCanvas(cDr);
                driftView.Plot(true);
                TLine* dLl = new TLine(-1.5, -16.0, -1.5, 16.0);
                TLine* dRl = new TLine( 1.5, -16.0,  1.5, 16.0);
                dLl->SetLineColor(kGray + 2); dLl->SetLineWidth(2); dLl->Draw();
                dRl->SetLineColor(kGray + 2); dRl->SetLineWidth(2); dRl->Draw();
                cDr->SaveAs((out_dir + "/diffusion_drift_lines.pdf").c_str());
                delete cDr;
            }
            delete hEnd;
        }
        fsum.close();

        // ---- comparison plot: measured vs tabulated D_T (and D_L) vs O2 % ----
        TCanvas* cCmp = new TCanvas("cCmp", "Diffusion vs O2", 900, 650);
        cCmp->SetGrid();
        TMultiGraph* mgD = new TMultiGraph();
        mgD->SetTitle("Diffusion vs O_{2} contamination (67/33 CF_{4}/iC_{4}H_{10});O_{2} [%];Diffusion coefficient [cm^{1/2}]");
        grDTmeas->SetMarkerStyle(20); grDTmeas->SetMarkerColor(kRed+1);  grDTmeas->SetLineColor(kRed+1);  grDTmeas->SetLineWidth(2); grDTmeas->SetTitle("D_{T} measured (drift fit)");
        grDTtab ->SetMarkerStyle(24); grDTtab ->SetMarkerColor(kRed+1);  grDTtab ->SetLineColor(kRed+1);  grDTtab ->SetLineStyle(2); grDTtab ->SetTitle("D_{T} tabulated (gas table)");
        grDLmeas->SetMarkerStyle(21); grDLmeas->SetMarkerColor(kBlue+1); grDLmeas->SetLineColor(kBlue+1); grDLmeas->SetLineWidth(2); grDLmeas->SetTitle("D_{L} measured (drift fit)");
        grDLtab ->SetMarkerStyle(25); grDLtab ->SetMarkerColor(kBlue+1); grDLtab ->SetLineColor(kBlue+1); grDLtab ->SetLineStyle(2); grDLtab ->SetTitle("D_{L} tabulated (gas table)");
        mgD->Add(grDTmeas,"LP"); mgD->Add(grDTtab,"LP");
        mgD->Add(grDLmeas,"LP"); mgD->Add(grDLtab,"LP");
        mgD->Draw("A");
        cCmp->BuildLegend(0.55, 0.15, 0.88, 0.40);
        cCmp->SaveAs((out_dir + "/diffusion_vs_o2.pdf").c_str());

        std::cout << "\n-> Wrote diffusion_vs_o2.txt, diffusion_vs_o2.pdf, "
                     "per-gas profiles and sigma plots to " << out_dir << "\n";
    } else if (switch_mode == 3) {
        // ==================================================================
        // MODE 3: SPATIAL SWEEP WITH READOUT RESPONSE  (single gas)   [v2]
        //
        // Fixes over v1:
        //   - sigma_x / sigma_t are now per-event (v1 pooled all x
        //     positions, so it measured the spawn spread, not diffusion)
        //   - N_EVT events per x position (v1 had one, so the adjusted
        //     curve was quantised in 1/21 = 4.76% steps)
        //   - thresholds expressed as a FRACTION of full cluster charge,
        //     so the scan does not depend on GAIN and N_PRIM
        //   - readout strip phase randomised per event
        //   - optional deposit (Landau) and gain (Polya) fluctuations,
        //     without which the adjusted curve must start at 100%
        // ==================================================================


        std::cout << "\n=========================================\n";
        std::cout << "MODE 3: SWEEP + READOUT RESPONSE  [v2]\n";
        std::cout << "=========================================\n";


        // ------------------- CONFIG -------------------
        const std::string gas_file_m3 = "gas_tables/cf4_20_iso_80.gas";


        const double READOUT_PITCH = 0.254;   // cm   <<CHECK>>
        const int    N_STRIPS      = 25;
        const double GAIN          = 1000.0;  //      <<CHECK>>
        const int    N_PRIM        = 30;      //      <<CHECK>> primaries per cluster


        const int    N_X           = 15;      // spawn x positions
        const int    N_EVT         = 10;      // events per x position
        const double X_MIN = -1.4, X_MAX = 1.4;


        // Default to the window the real data covers (0-11 cm).
        // Set FULL_RANGE for the whole cage.
        const bool   FULL_RANGE = false;
        const double Y_MIN  = FULL_RANGE ? -15.0 : -15.5;
        const double Y_MAX  = FULL_RANGE ?  15.0 :  -5.0;
        const double y_step = 0.1;


        const double Y_ACCEPT = -15.98;
        const double X_ACCEPT =   1.50;


        // Thresholds as a fraction of full cluster charge (N_PRIM * GAIN).
        // v1 showed the interesting region is 0.27 - 0.53, so sample it well.
        const std::vector<double> THR_FRAC =
            {0.15, 0.25, 0.30, 0.35, 0.40, 0.45, 0.55};
        const double THR_NOMINAL_FRAC = 0.35;


        // Fluctuations.  With both off, every cluster is identical and the
        // adjusted curve is forced to start at 100%; the real curve starts
        // at 90%, so leave these on for any comparison to data.
        const bool   FLUCT_DEPOSIT = true;
        const double LANDAU_WIDTH  = 0.15 * N_PRIM;   // <<CHECK>>
        const bool   FLUCT_GAIN    = true;
        const double POLYA_F       = 0.5;  // rel. variance of single-e gain <<CHECK>>


        // Ramo / induced-signal machinery removed: we only count electrons
        // reaching the M-THGEM. The weighting field it used came from an
        // idealized ComponentAnalyticField (infinite planes at +/-15), not
        // from our Elmer solve, so it was the one place a field was invented
        // outside Gmsh/Elmer. Efficiency does not depend on it.


        if (!gas.LoadGasFile(gas_file_m3)) {
            std::cerr << "FATAL: could not load " << gas_file_m3 << "\n";
            gSystem->Exit(1);
        }
        std::cout << "-> Gas: " << gas_file_m3 << "\n";


        const int nY = int(std::lround((Y_MAX - Y_MIN) / y_step)) + 1;
        const int evtPerY = N_X * N_EVT;
        std::cout << "-> " << nY << " heights x " << N_X << " positions x "
                  << N_EVT << " events x " << N_PRIM << " primaries\n";
        std::cout << "-> " << evtPerY << " events per height (resolution "
                  << 100.0 / evtPerY << "%)\n";
        std::cout << "-> " << long(nY) * evtPerY * N_PRIM << " drift lines total\n\n";


        const double Q_FULL = N_PRIM * GAIN;
        std::vector<double> THR_E;
        for (double f : THR_FRAC) THR_E.push_back(f * Q_FULL);
        const double THR_NOM = THR_NOMINAL_FRAC * Q_FULL;


        AvalancheMC m3drift;
        m3drift.SetSensor(&sensor);
        m3drift.SetDistanceSteps(0.02);
        m3drift.SetTimeWindow(0.0, 500000.0);
        m3drift.EnableAttachment();
        m3drift.EnableDiffusion();     // ← ADD

        std::ofstream fMain(out_dir + "/mode3_curves.txt");
        fMain << "Spawn_Y\tDrift_L\tUnadj_pct\tAdj_pct\tsigma_x\tsigma_t"
                 "\tmean_peak_frac\tn_events\n";
        std::ofstream fScan(out_dir + "/mode3_threshold_scan.txt");
        fScan << "Spawn_Y\tDrift_L";
        for (double f : THR_FRAC) fScan << "\tfrac_" << f;
        fScan << "\n";


        TGraph* grUnadj = new TGraph();
        TGraph* grAdj   = new TGraph();
        TGraph* grSigX  = new TGraph();
        TGraph* grSigT  = new TGraph();
        TGraph* grFrac  = new TGraph();
        std::vector<TGraph*> grScan;
        for (size_t k = 0; k < THR_FRAC.size(); ++k) grScan.push_back(new TGraph());


        int pt = 0;
        for (int iy = 0; iy < nY; ++iy) {
            const double y_spawn = Y_MIN + y_step * iy;
            const double drift_L = y_spawn + 15.0;


            long n_launched = 0, n_arrived = 0;
            int  n_events = 0, n_fired = 0;
            std::vector<int> fired_scan(THR_FRAC.size(), 0);
            double sigx_sum = 0, sigt_sum = 0, frac_sum = 0;
            int    sig_n = 0;


            for (int ix = 0; ix < N_X; ++ix) {
                const double x0 = (N_X == 1) ? 0.0
                    : X_MIN + (X_MAX - X_MIN) * ix / double(N_X - 1);


                for (int iev = 0; iev < N_EVT; ++iev) {


                    int n_this = N_PRIM;
                    if (FLUCT_DEPOSIT) {
                        n_this = int(std::lround(
                                     gRandom->Landau(N_PRIM, LANDAU_WIDTH)));
                        if (n_this < 1)          n_this = 1;
                        if (n_this > 10 * N_PRIM) n_this = 10 * N_PRIM;
                    }


                    // random strip phase: real clusters land anywhere
                    // relative to the strip pattern
                    const double phase =
                        (RndmUniform() - 0.5) * READOUT_PITCH;


                    std::vector<int> strip_count(N_STRIPS, 0);
                    std::vector<double> ax, at;
                    ax.reserve(n_this); at.reserve(n_this);


                    for (int e = 0; e < n_this; ++e) {
                        ++n_launched;
                        m3drift.DriftElectron(x0, y_spawn, 0.005, 0.0);
                        const int nPoints = m3drift.GetNumberOfDriftLinePoints();
                        if (nPoints <= 0) continue;
                        double x1 = 0, y1 = 0, z1 = 0, t1 = 0;
                        m3drift.GetDriftLinePoint(nPoints - 1, x1, y1, z1, t1);
                        if (y1 < Y_ACCEPT && std::abs(x1) < X_ACCEPT) {
                            ++n_arrived;
                            ax.push_back(x1);
                            at.push_back(t1);
                            int is = int(std::lround((x1 + phase) / READOUT_PITCH))
                                     + N_STRIPS / 2;
                            if (is >= 0 && is < N_STRIPS) strip_count[is]++;
                        }
                    }


                    // ---- per-event spreads (this is the diffusion sigma) ----
                    if (ax.size() > 2) {
                        double mx = 0, mt = 0;
                        for (size_t j = 0; j < ax.size(); ++j) { mx += ax[j]; mt += at[j]; }
                        mx /= ax.size(); mt /= at.size();
                        double vx = 0, vt = 0;
                        for (size_t j = 0; j < ax.size(); ++j) {
                            vx += (ax[j] - mx) * (ax[j] - mx);
                            vt += (at[j] - mt) * (at[j] - mt);
                        }
                        sigx_sum += std::sqrt(vx / (ax.size() - 1));
                        sigt_sum += std::sqrt(vt / (at.size() - 1));
                        ++sig_n;
                    }


                    // ---- readout response ----
                    int peak = 0;
                    for (int s = 0; s < N_STRIPS; ++s)
                        peak = std::max(peak, strip_count[s]);


                    double q_peak = peak * GAIN;
                    if (FLUCT_GAIN && peak > 0) {
                        q_peak += gRandom->Gaus(0.0,
                                     GAIN * std::sqrt(POLYA_F * peak));
                        if (q_peak < 0.0) q_peak = 0.0;
                    }


                    ++n_events;
                    if (q_peak > THR_NOM) ++n_fired;
                    for (size_t k = 0; k < THR_E.size(); ++k)
                        if (q_peak > THR_E[k]) fired_scan[k]++;
                    if (!ax.empty()) frac_sum += double(peak) / ax.size();
                }
            }


            const double unadj = 100.0 * n_arrived / std::max(1L, n_launched);
            const double adj   = 100.0 * n_fired   / std::max(1,  n_events);
            const double sigx  = sig_n ? sigx_sum / sig_n : 0.0;
            const double sigt  = sig_n ? sigt_sum / sig_n : 0.0;
            const double mfrac = frac_sum / std::max(1, n_events);


            fMain << y_spawn << "\t" << drift_L << "\t" << unadj << "\t" << adj
                  << "\t" << sigx << "\t" << sigt << "\t" << mfrac
                  << "\t" << n_events << "\n";
            fScan << y_spawn << "\t" << drift_L;
            for (size_t k = 0; k < THR_E.size(); ++k)
                fScan << "\t" << 100.0 * fired_scan[k] / std::max(1, n_events);
            fScan << "\n";


            grUnadj->SetPoint(pt, drift_L, unadj);
            grAdj  ->SetPoint(pt, drift_L, adj);
            grSigX ->SetPoint(pt, drift_L, sigx);
            grSigT ->SetPoint(pt, drift_L, sigt);
            grFrac ->SetPoint(pt, drift_L, 100.0 * mfrac);
            for (size_t k = 0; k < THR_E.size(); ++k)
                grScan[k]->SetPoint(pt, drift_L,
                                    100.0 * fired_scan[k] / std::max(1, n_events));
            ++pt;


            std::cout << "  L = " << std::setw(5) << std::fixed
                      << std::setprecision(1) << drift_L
                      << "  unadj " << std::setw(5) << unadj
                      << "%  adj " << std::setw(5) << adj
                      << "%  sx " << std::setprecision(4) << sigx
                      << "  st " << std::setprecision(1) << sigt << " ns\n"
                      << std::flush;
        }
        fMain.close();
        fScan.close();


        // ------------------- plots -------------------
        gStyle->SetOptStat(0);
        int cols[] = {kBlack, kBlue, kGreen + 2, kOrange + 1, kRed,
                      kMagenta + 2, kAzure + 2};


        TCanvas* c1 = new TCanvas("c_m3", "readout response", 1200, 800);
        TMultiGraph* mg = new TMultiGraph();
        grUnadj->SetLineColor(kBlack); grUnadj->SetLineWidth(3);
        grUnadj->SetTitle("Unadjusted (electrons arriving)");
        grAdj->SetLineColor(kRed + 1); grAdj->SetLineWidth(3);
        grAdj->SetTitle(Form("Adjusted (threshold %.0f%% of cluster)",
                             100 * THR_NOMINAL_FRAC));
        mg->Add(grUnadj, "L"); mg->Add(grAdj, "L");
        mg->SetTitle("Pass fraction vs drift distance;"
                     "Distance from M-THGEM [cm];Passed [%]");
        mg->Draw("AL");
        mg->GetYaxis()->SetRangeUser(0, 105);
        c1->BuildLegend(0.50, 0.74, 0.88, 0.88);
        c1->SaveAs((out_dir + "/mode3_pass_vs_L.pdf").c_str());


        TCanvas* c2 = new TCanvas("c_m3s", "threshold scan", 1200, 800);
        TMultiGraph* mgs = new TMultiGraph();
        for (size_t k = 0; k < THR_FRAC.size(); ++k) {
            grScan[k]->SetLineColor(cols[k % 7]);
            grScan[k]->SetLineWidth(2);
            grScan[k]->SetTitle(Form("threshold %.0f%% of cluster",
                                     100 * THR_FRAC[k]));
            mgs->Add(grScan[k], "L");
        }
        mgs->SetTitle("Efficiency vs drift distance, threshold scan;"
                      "Distance from M-THGEM [cm];Events over threshold [%]");
        mgs->Draw("AL");
        mgs->GetYaxis()->SetRangeUser(0, 105);
        c2->BuildLegend(0.55, 0.60, 0.88, 0.88);
        c2->SaveAs((out_dir + "/mode3_threshold_scan.pdf").c_str());


        TCanvas* c3 = new TCanvas("c_m3d", "spreads", 1200, 500);
        c3->Divide(2, 1);
        c3->cd(1);
        grSigX->SetLineColor(kBlue); grSigX->SetLineWidth(3);
        grSigX->SetMarkerStyle(20);
        grSigX->SetTitle("Per-cluster transverse spread;"
                         "Distance from M-THGEM [cm];#sigma_{x} [cm]");
        grSigX->Draw("ALP");
        c3->cd(2);
        grSigT->SetLineColor(kRed + 1); grSigT->SetLineWidth(3);
        grSigT->SetMarkerStyle(20);
        grSigT->SetTitle("Per-cluster arrival time spread;"
                         "Distance from M-THGEM [cm];#sigma_{t} [ns]");
        grSigT->Draw("ALP");
        c3->SaveAs((out_dir + "/mode3_spreads.pdf").c_str());


        TCanvas* c4 = new TCanvas("c_m3f", "peak fraction", 1200, 800);
        grFrac->SetLineColor(kMagenta + 2); grFrac->SetLineWidth(3);
        grFrac->SetTitle("Charge on peak strip;"
                         "Distance from M-THGEM [cm];Peak strip / arrived [%]");
        grFrac->Draw("AL");
        c4->SaveAs((out_dir + "/mode3_peak_fraction.pdf").c_str());


        std::cout << "\n-> Mode 3 v2 outputs written to " << out_dir << "\n";
    } else if (switch_mode == 4) {
        // ==================================================================
        // MODE 4 - RIPPLE TEST: does measured D_T depend on x0?
        //
        // Mode 2 found measured D_T ~1.3% above the Magboltz tabulated value.
        // Hypothesis: the field ripple from the periodic strips at x=+/-1.5
        // adds a transverse kick on top of thermal diffusion, and the ripple
        // is strongest near the strips, weakest on the center axis.
        //
        // Test: measure D_T from a cloud launched at several FIXED x0 columns,
        // from x0=0 (deep center, minimal ripple) out toward x0=1.3 (next to
        // the strips, max ripple). Each column keeps its full statistics (the
        // launch x0 does not change during the fit because we subtract the
        // per-plane mean before taking the variance -- sigma_x measures spread
        // about the local centroid, not absolute position).
        //
        // If D_T(x0) rises toward the strips, the ripple is proven. If it is
        // flat across x0, the excess is something else (fit/step systematic).
        // Attachment OFF: pure diffusion only.
        // ==================================================================
        const std::string gas_file_m4 = "gas_tables_v2/cf4_67_iso_33_o2_0p000pct.gas";
        if (!gas.LoadGasFile(gas_file_m4)) {
            std::cerr << "FATAL: could not load " << gas_file_m4 << "\n";
            gSystem->Exit(1);
        }
        std::cout << "-> Gas: " << gas_file_m4 << "\n";

        const int    nElectrons = 20000;   // per x0 column; high stats for a clean per-column D_T
        const double y_start     =  0.0;
        const double z_start     =  0.005;
        const double y_stop      = -15.0;
        const double dy_check    =  0.5;

        // x0 columns to probe, from center to just inside the strip plane.
        std::vector<double> x0_cols = {0.0, 0.3, 0.6, 0.9, 1.1, 1.3};

        // Field (gas-independent), same window as mode 2 so the quoted E matches.
        double sum_e4 = 0.0; int n_e4 = 0;
        for (double x = -1.490; x <= 1.491; x += 0.1)
            for (double y = -15.0; y <= 15.0; y += 0.2) {
                double ex=0,ey=0,ez=0,v=0; Medium* m=nullptr; int st=0;
                elm.ElectricField(x,y,0.005,ex,ey,ez,v,m,st);
                if (st==0){ sum_e4 += std::sqrt(ex*ex+ey*ey+ez*ez); ++n_e4; }
            }
        const double avg_e = (n_e4>0)? sum_e4/n_e4 : 33.0;
        std::cout << "-> Average field for mode 4: " << avg_e << " V/cm\n";

        // Tabulated D_T for this gas/field (the reference every column is compared to).
        double dl_tab=0, dt_tab=0;
        gas.ElectronDiffusion(0.0, avg_e, 0.0, 0.,0.,0., dl_tab, dt_tab);
        std::cout << "-> Tabulated D_T (reference) = " << dt_tab << " cm^1/2\n";

        const int nCheck = int(std::lround((y_start - y_stop)/dy_check));
        std::vector<double> yCheck(nCheck);
        for (int c=0;c<nCheck;++c) yCheck[c] = y_start - dy_check*(c+1);

        std::ofstream fr(out_dir + "/ripple_test_DT_vs_x0.txt");
        fr << "x0_cm\tDT_meas\tDT_err\tDT_tab\tratio\tmean_abs_x_during_drift\n";

        TGraph* grDTx = new TGraph();   // measured D_T vs x0

        for (int ix=0; ix<(int)x0_cols.size(); ++ix) {
            const double x0 = x0_cols[ix];

            std::vector<long long> nAt(nCheck,0);
            std::vector<double> sX(nCheck,0.0), sX2(nCheck,0.0);
            double sumAbsX = 0.0; long nAbsX = 0;    // track where the cloud actually sits

            AvalancheMC dif;
            dif.SetSensor(&sensor);
            dif.SetDistanceSteps(0.01);
            dif.SetTimeWindow(0.0, 500000.0);
            dif.EnableDiffusion();
            // attachment OFF -- pure diffusion

            std::cout << "\n  x0 = " << x0 << " cm : drifting " << nElectrons << " ...\n";
            for (int e=0;e<nElectrons;++e) {
                dif.DriftElectron(x0, y_start, z_start, 0.0);
                const int nPts = dif.GetNumberOfDriftLinePoints();
                if (nPts < 2) continue;
                int c=0; double xP,yP,zP,tP;
                dif.GetDriftLinePoint(0,xP,yP,zP,tP);
                for (int p=1;p<nPts && c<nCheck;++p) {
                    double x,y,z,t; dif.GetDriftLinePoint(p,x,y,z,t);
                    while (c<nCheck && yP>yCheck[c] && y<=yCheck[c]) {
                        const double dyseg = yP - y;
                        const double f = (dyseg>1e-12)? (yP-yCheck[c])/dyseg : 0.0;
                        const double xi = xP + f*(x-xP);
                        nAt[c]+=1; sX[c]+=xi; sX2[c]+=xi*xi;
                        sumAbsX += std::abs(xi); ++nAbsX;
                        ++c;
                    }
                    xP=x; yP=y; zP=z; tP=t;
                }
                if ((e+1)%2000==0) std::cout << "    " << (e+1) << " / " << nElectrons << "\r" << std::flush;
            }
            std::cout << "\n";

            // sigma_x about the per-plane centroid (so absolute x0 does not inflate it),
            // fit sigma_x = D_T * sqrt(L) through the origin.
            TGraph* grSig = new TGraph(); int k=0;
            for (int c=0;c<nCheck;++c) {
                if (nAt[c]<50) continue;
                const double N=double(nAt[c]);
                const double L=y_start-yCheck[c];
                const double mx=sX[c]/N;
                double vx=sX2[c]/N-mx*mx; if (vx<0) vx=0;
                grSig->SetPoint(k++, std::sqrt(L), std::sqrt(vx));
            }
            TF1* fT = new TF1(Form("fT4_%d",ix), "[0]*x", 0.0, std::sqrt(y_start-y_stop)*1.05);
            fT->SetParameter(0, dt_tab); grSig->Fit(fT,"Q");
            const double DT_fit=fT->GetParameter(0), DT_err=fT->GetParError(0);
            const double meanAbsX = (nAbsX>0)? sumAbsX/nAbsX : x0;

            std::cout << "    x0=" << x0 << "  measured D_T=" << DT_fit << " +/- " << DT_err
                      << "  ratio=" << DT_fit/dt_tab << "  <|x|>=" << meanAbsX << " cm\n";
            fr << x0 << "\t" << DT_fit << "\t" << DT_err << "\t" << dt_tab << "\t"
               << DT_fit/dt_tab << "\t" << meanAbsX << "\n";
            grDTx->SetPoint(ix, x0, DT_fit);
            delete grSig;
        }
        fr.close();

        // Plot: measured D_T vs x0, with the tabulated reference as a dashed line.
        TCanvas* cR = new TCanvas("cR","Ripple test",900,650);
        cR->SetGrid();
        grDTx->SetTitle("Ripple test: measured D_{T} vs launch position x_{0};x_{0} [cm];Measured D_{T} [cm^{1/2}]");
        grDTx->SetMarkerStyle(20); grDTx->SetMarkerSize(1.1);
        grDTx->SetMarkerColor(kRed+1); grDTx->SetLineColor(kRed+1); grDTx->SetLineWidth(2);
        grDTx->Draw("ALP");
        TLine* lt = new TLine(0.0, dt_tab, x0_cols.back(), dt_tab);
        lt->SetLineColor(kBlack); lt->SetLineStyle(2); lt->SetLineWidth(2); lt->Draw();
        TLegend* lg = new TLegend(0.15,0.75,0.6,0.88); lg->SetBorderSize(0);
        lg->AddEntry(grDTx,"measured D_{T}(x_{0})","lp");
        lg->AddEntry(lt,Form("tabulated D_{T} = %.4f",dt_tab),"l");
        lg->Draw();
        cR->SaveAs((out_dir + "/ripple_test_DT_vs_x0.pdf").c_str());

        std::cout << "\n-> Mode 4 ripple test written to " << out_dir
                  << " (ripple_test_DT_vs_x0.txt and .pdf)\n";
        std::cout << "   If D_T rises toward x0=1.3 (the strips), the ripple is confirmed.\n";
    } else if (switch_mode == 6) {
        // ==================================================================
        // MODE 6 - FULL E-FIELD DUMP TO CSV
        //
        // Exact electric field on a dense grid over the whole domain:
        //   x in [-1.6, 1.6], y in [-16, 16], at z = 0.005 cm (the drift plane).
        //
        // Probes elm.ElectricField DIRECTLY (not the Sensor), because the
        // Sensor's SetArea clips to |x|<=1.5 and we want values OUTSIDE the
        // electrode planes too (out to +/-1.6).
        //
        // Every grid point is written, including points inside metal / off-mesh
        // (flagged by status != 0), so the CSV is complete and filterable.
        // Columns: x, y, z, ex, ey, ez, emag, v, status
        //   ex = transverse (x) field component  [V/cm]
        //   ey = drift (y) field component        [V/cm]
        //   status = 0 -> valid gas point; != 0 -> inside metal or off-mesh
        // ==================================================================
        const double X_LO = -1.6, X_HI = 1.6;
        const double Y_LO = -16.0, Y_HI = 16.0;
        const double Z_PROBE = 0.005;
        const double DX = 0.01;   // 100 um grid in x
        const double DY = 0.01;   // 100 um grid in y

        const int nx = int(std::lround((X_HI - X_LO) / DX)) + 1;
        const int ny = int(std::lround((Y_HI - Y_LO) / DY)) + 1;
        std::cout << "\n=========================================\n";
        std::cout << "MODE 6: FULL E-FIELD DUMP\n";
        std::cout << "  x [" << X_LO << ", " << X_HI << "] step " << DX << "  (" << nx << " cols)\n";
        std::cout << "  y [" << Y_LO << ", " << Y_HI << "] step " << DY << "  (" << ny << " rows)\n";
        std::cout << "  z = " << Z_PROBE << " cm\n";
        std::cout << "  total points = " << (long)nx * ny << "\n";
        std::cout << "=========================================\n";

        std::ofstream fcsv(out_dir + "/efield_full_strip.csv");
        fcsv << "x,y,z,ex,ey,ez,emag,v,status\n";
        fcsv << std::setprecision(10);

        long written = 0, valid = 0;
        for (int ix = 0; ix < nx; ++ix) {
            const double x = X_LO + DX * ix;
            for (int iy = 0; iy < ny; ++iy) {
                const double y = Y_LO + DY * iy;
                double ex = 0., ey = 0., ez = 0., v = 0.;
                Medium* m = nullptr;
                int status = 0;
                elm.ElectricField(x, y, Z_PROBE, ex, ey, ez, v, m, status);
                const double emag = std::sqrt(ex*ex + ey*ey + ez*ez);
                fcsv << x << "," << y << "," << Z_PROBE << ","
                     << ex << "," << ey << "," << ez << ","
                     << emag << "," << v << "," << status << "\n";
                ++written;
                if (status == 0) ++valid;
            }
            if ((ix + 1) % 20 == 0)
                std::cout << "  column " << (ix + 1) << " / " << nx << "\r" << std::flush;
        }
        fcsv.close();
        std::cout << "\n-> Wrote " << written << " rows (" << valid << " valid gas points) to "
                  << out_dir << "/efield_full_strip.csv\n";
    } else if (switch_mode == 7) {
        // ==================================================================
        // MODE 7 - MICRO ELECTRIC FIELD PROBE  (3x3 "crafting table")
        //
        // Nine 1 cm x 1 cm squares tiled over the detector so we can see the
        // transverse field ripple with a continuous, high-resolution map.
        // Layout (columns = x, rows = y), tiling the full 3 cm width at three
        // heights so both electrode planes AND the uniform center are sampled:
        //
        //        x: [+0.5,+1.5]   [-0.5,+0.5]   [-1.5,-0.5]
        //           (right plane) (  center  )  (left plane)
        //   y top    r2c0           r2c1          r2c2        [+12,+13]
        //   y mid    r1c0           r1c1          r1c2        [-0.5,+0.5]
        //   y bot    r0c0           r0c1          r0c2        [-13,-12]
        //
        // Each square is sampled on an N_PER_SIDE x N_PER_SIDE grid.
        //   N_PER_SIDE = 2000 -> 5 um spacing, 4.0 M points/square, 36 M total.
        // The Elmer mesh is ~200 um near the electrodes, so 5 um oversamples the
        // field enough to look perfectly continuous without inventing physics.
        // Bump N_PER_SIDE if you want finer (e.g. 10000 -> 1 um, 1e8/square),
        // but note the CSVs grow ~linearly with the point count.
        //
        // Probes elm.ElectricField DIRECTLY (not the Sensor) so values right at
        // and behind the electrode planes (|x| up to 1.5) are included.
        // One CSV per square: efield_micro_strip_c{C}.csv
        // Columns: x, y, z, ex, ey, ez, emag, v, status
        //   ex = transverse (x) field  [V/cm];  ey = drift (y) field  [V/cm]
        //   status = 0 -> valid gas point;  != 0 -> inside metal or off-mesh
        // ==================================================================
        // REDUCED TO 3 SQUARES (one row across the width) to save disk.
        // All at the clean middle y-band; they sweep x so we capture:
        //   c0 right electrode plane | c1 uniform center | c2 left electrode plane
        // The three y-heights were only a uniformity check already confirmed by
        // the MODE 6 full dump, so one height loses nothing important.
        const int    N_PER_SIDE = 500;       // 500 -> ~20 um grid, 250k pts/sq, ~25 MB CSV (tunable; 2000 -> 5 um)
        const double SQ         = 1.0;       // square side [cm]
        const double Z_PROBE    = 0.005;     // drift plane
        const double step       = SQ / (N_PER_SIDE - 1);
        const double Y0         = -0.5;      // single y-band: clean middle, y in [-0.5, 0.5]

        // 3 square lower-left x corners: right plane, center, left plane
        const double x_lo[3] = { 0.5, -0.5, -1.5 };

        std::cout << "\n=========================================\n";
        std::cout << "MODE 7: MICRO E-FIELD PROBE (3 squares across the width)\n";
        std::cout << "  3 squares, each " << SQ << " x " << SQ << " cm, y in ["
                  << Y0 << "," << Y0+SQ << "]\n";
        std::cout << "  grid " << N_PER_SIDE << " x " << N_PER_SIDE
                  << "  (" << step*1e4 << " um spacing, "
                  << (long)N_PER_SIDE*N_PER_SIDE << " pts/square)\n";
        std::cout << "  total points = " << 3L*N_PER_SIDE*N_PER_SIDE << "\n";
        std::cout << "  z = " << Z_PROBE << " cm\n";
        std::cout << "=========================================\n";

        long grand_total = 0, grand_valid = 0;
        for (int col = 0; col < 3; ++col) {
            const double X0 = x_lo[col];
            std::ostringstream fn;
            fn << out_dir << "/efield_micro_strip_c" << col << ".csv";
            std::ofstream fcsv(fn.str());
            fcsv << "x,y,z,ex,ey,ez,emag,v,status\n";
            fcsv << std::setprecision(10);

            long written = 0, valid = 0;
            for (int ix = 0; ix < N_PER_SIDE; ++ix) {
                const double x = X0 + step * ix;
                for (int iy = 0; iy < N_PER_SIDE; ++iy) {
                    const double y = Y0 + step * iy;
                    double ex = 0., ey = 0., ez = 0., v = 0.;
                    Medium* m = nullptr;
                    int status = 0;
                    elm.ElectricField(x, y, Z_PROBE, ex, ey, ez, v, m, status);
                    const double emag = std::sqrt(ex*ex + ey*ey + ez*ez);
                    fcsv << x << "," << y << "," << Z_PROBE << ","
                         << ex << "," << ey << "," << ez << ","
                         << emag << "," << v << "," << status << "\n";
                    ++written;
                    if (status == 0) ++valid;
                }
                if ((ix + 1) % 100 == 0)
                    std::cout << "  c" << col
                              << "  col " << (ix + 1) << " / " << N_PER_SIDE << "\r" << std::flush;
            }
            fcsv.close();
            grand_total += written; grand_valid += valid;
            std::cout << "\n  -> c" << col << "  x[" << X0 << "," << X0+SQ
                      << "] y[" << Y0 << "," << Y0+SQ << "]  "
                      << written << " rows (" << valid << " valid) -> " << fn.str() << "\n";
        }
        std::cout << "\n== MODE 7 done: " << grand_total << " points total ("
                  << grand_valid << " valid gas points) across 3 squares ==\n";
    }




    gSystem->Exit(0);
    return 0;
}