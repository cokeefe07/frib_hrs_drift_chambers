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
#include <TH2D.h>
#include <TEllipse.h>
#include <TBox.h>
#include <TStyle.h>
#include <TPad.h>


// OpenMP intentionally not used: parallel AvalancheMC drifts caused
// cross-thread corruption. All modes run serial (matches strip_cage.cc).


using namespace Garfield;




// =====================================================================
//  E-FIELD POINT MAP  --  grid of coloured dots, diameter = (3/4)*spacing
//  Ported from strip_cage.cc. Reference lines mark the WIRE planes.
//
//  ASPECT: circles have radius 0.375*d in DATA units, so the canvas aspect
//  is derived from the window you request. Zoomed windows read best.
//
//  zmin/zmax clip the colour range (0,0 = autoscale). Autoscale is dominated
//  by the ~16x field concentration right at the wires, which flattens the
//  whole drift region to one colour -- clip it to see the bulk structure.
// =====================================================================
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


  std::cout << "  field point map (" << tag << "): " << px.size()
            << " points, |E| range " << emin << " .. " << emax << " V/cm\n";


  gStyle->SetPalette(kBird);
  const int nCol = TColor::GetNumberOfColors();


  int canvasH = int(canvasW * (yhi - ylo) / (xhi - xlo));
  if (canvasH < 300)  canvasH = 300;
  if (canvasH > 2400) canvasH = 2400;


  TCanvas* c = new TCanvas(("cFieldPts_" + tag).c_str(),
                           "E field point map", canvasW, canvasH);
  c->SetRightMargin(0.18);
  c->SetLeftMargin(0.14);


  TH2D* hF = new TH2D(("hFieldFrame_" + tag).c_str(),
                      ";x [cm];y [cm]", 10, xlo, xhi, 10, ylo, yhi);
  hF->SetStats(0);
  hF->SetBinContent(1, 1, emin);      // one entry so ROOT draws the palette
  hF->SetMinimum(emin);
  hF->SetMaximum(emax);
  hF->GetZaxis()->SetTitle("|E| [V/cm]");
  hF->Draw("COLZ");
  gPad->Update();


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


  // wire planes for reference
  TLine* wL = new TLine(-1.5, ylo, -1.5, yhi);
  TLine* wR = new TLine( 1.5, ylo,  1.5, yhi);
  wL->SetLineColor(kGray + 2); wL->SetLineWidth(2); wL->Draw();
  wR->SetLineColor(kGray + 2); wR->SetLineWidth(2); wR->Draw();


  gPad->RedrawAxis();
  c->SaveAs((out_dir + "/field_point_map_" + tag + ".pdf").c_str());
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
    gas.LoadGasFile("gas_tables/cf4_67_iso_33.gas");


    ComponentElmer elm("field_cage/mesh.header", "field_cage/mesh.elements", "field_cage/mesh.nodes",
                       "dielectrics.dat", "field_cage/field_cage.result", "cm");
    elm.SetMedium(0, &gas);
    elm.EnablePeriodicityZ();


    Sensor sensor;
    sensor.AddComponent(&elm);
    sensor.SetArea(-1.5, -16.0, -100.0, 1.5, 16.0, 100);


    AvalancheMC drift;
    drift.SetSensor(&sensor);
    drift.SetDistanceSteps(0.02);


    drift.SetTimeWindow(0.0, 500000.0);
    drift.EnableAttachment();
    drift.EnableDiffusion();


    ViewDrift driftView;
    driftView.SetArea(-2.0, -16.0, -1.0, 2.0, 16.0, 1.0);




    // ======================================================================
    // SWITCH: 0 = Visual Diagnostic (67/33), 1 = Multi-Gas Spatial Sweep,
    //         6 = Full E-field dump,  7 = Micro E-field probe (3x3 grid)
    // ======================================================================
    int switch_mode = 0;


    if (switch_mode == 0) {
        std::cout << "\n=========================================\n";
        std::cout << "MODE 0: VISUAL DIAGNOSTIC (67% CF4)\n";
        std::cout << "=========================================\n";


        gas.LoadGasFile("gas_tables/cf4_67_iso_33.gas");
        drift.EnablePlotting(&driftView);


        // ---- deterministic grid seeding: n electrons per layer, stepped in y ----
        const int    n_per_layer = 10;      // electrons spread across x at each layer
        const double y_step_vis  = 1;    // cm between layers  (31 cm / 2.0 -> 16 layers)
        const double y_first     = -15.0;
        const double y_last      =  15.0;
        const double x_lo        = -1.4;   // same span MODE 1 uses
        const double x_hi        =  1.4;


        int nTotal = 0, nCollected = 0;
        std::cout << "Drifting " << n_per_layer << " electrons per layer, y from "
                  << y_first << " to " << y_last
                  << " cm in bsteps of " << y_step_vis << " cm\n";


        for (double y0 = y_first; y0 <= y_last + 1e-9; y0 += y_step_vis) {
            int layerCollected = 0;
            for (int e = 0; e < n_per_layer; ++e) {
                const double x0 = (n_per_layer == 1)
                                ? 0.5 * (x_lo + x_hi)
                                : x_lo + (x_hi - x_lo) * e / (n_per_layer - 1.0);
                const double z0 = 0.005, t0 = 0.0;


                auto te0 = std::chrono::steady_clock::now();
                drift.DriftElectron(x0, y0, z0, t0);
                auto te1 = std::chrono::steady_clock::now();
                ++nTotal;


                const int nPoints = drift.GetNumberOfDriftLinePoints();
                double x1 = 0, y1 = 0, z1 = 0, t1 = 0;
                if (nPoints > 0) drift.GetDriftLinePoint(nPoints - 1, x1, y1, z1, t1);
                // Same acceptance as MODE 1: reached the M-THGEM at y=-16.
                if (nPoints > 0 && (y1 < -15.98) && (std::abs(x1) < 1.5)) {
                    ++layerCollected;
                    ++nCollected;
                }
                std::cout << "  x0=" << x0 << " y0=" << y0
                << "  npts=" << nPoints
                << "  end=(" << x1 << "," << y1 << ")"
                << "  t=" << t1 << " ns"
                << "  wall=" << std::chrono::duration<double>(te1 - te0).count() << " s\n"
                << std::flush;


            }
            std::cout << "  y = " << y0 << " cm : "
                      << layerCollected << "/" << n_per_layer
                      << " reached the M-THGEM\n";
        }


        std::cout << "RESULTS: " << nCollected << " / " << nTotal
                  << " reached the M-THGEM ("
                  << (100.0 * nCollected / nTotal) << "%)\n\n";


        // 1. Drift Lines
        TCanvas* c1 = new TCanvas("c1", "Electron Drift Lines", 600, 600);
        driftView.SetCanvas(c1);
        driftView.Plot(true);
        TLine* stripL = new TLine(-1.5, -16.0, -1.5, 16.0);
        TLine* stripR = new TLine( 1.5, -16.0,  1.5, 16.0);
        stripL->SetLineColor(kGray+2); stripL->SetLineWidth(2); stripL->Draw();
        stripR->SetLineColor(kGray+2); stripR->SetLineWidth(2); stripR->Draw();
        c1->SaveAs((out_dir + "/field_drift_lines_67_33.pdf").c_str());
        std::cout << "-> Saved " << out_dir << "/field_drift_lines_67_33.pdf\n";


        // 2. Field Map
        TCanvas* cMap = new TCanvas("cMap", "Electric Field Contours", 800, 800);
        ViewField fieldView;
        fieldView.SetSensor(&sensor);
        fieldView.SetCanvas(cMap);
        fieldView.SetPlane(0, 0, 1, 0, 0, 0.005);
        fieldView.SetArea(1.4, 14.759, 1.6, 14.959);
        fieldView.SetNumberOfContours(50);
        fieldView.SetNumberOfSamples2d(400, 400);
        fieldView.SetVoltageRange(-1002., -980.);
        fieldView.PlotContour("v");
        TLine* stripL2 = new TLine(-1.5, -16.0, -1.5, 16.0);
        TLine* stripR2 = new TLine( 1.5, -16.0,  1.5, 16.0);
        stripL2->SetLineColor(kGray+2); stripL2->SetLineWidth(2); stripL2->Draw();
        stripR2->SetLineColor(kGray+2); stripR2->SetLineWidth(2); stripR2->Draw();
        cMap->SaveAs((out_dir + "/field_map_contours_67_33.pdf").c_str());
        std::cout << "-> Saved " << out_dir << "/field_map_contours_67_33.pdf\n";


        // --- E-field point maps ---
        // Bulk drift region. Tight colour clip: the wire cage is uniform to ~1 V/cm,
        // so a wide range would render the whole thing a single colour.
        PlotFieldPointMap(sensor, out_dir, "zoom", 0.05,
                          -1.9, 1.9, -2.0, 2.0, 30.0, 38.0);


        // Top of the cage: wires stop at y = +14.859, and the leakage field the
        // wires must terminate peaks here.
        PlotFieldPointMap(sensor, out_dir, "topend", 0.05,
                          -1.9, 1.9, 12.0, 16.0, 25.0, 50.0);


        // Full height, coarse grid, tall canvas.
        PlotFieldPointMap(sensor, out_dir, "full", 0.25,
                          -1.9, 1.9, -15.9, 15.9, 30.0, 38.0, 0.005, 320);




        // 3. Drift Velocity Graph
        ViewMedium mediumView;
        mediumView.SetMedium(&gas);
        TCanvas* cVel = new TCanvas("cVel", "Drift Velocity", 800, 600);
        mediumView.SetCanvas(cVel);
        mediumView.PlotElectronVelocity('e');
        cVel->SaveAs((out_dir + "/field_drift_velocity_67_33.pdf").c_str());
        std::cout << "-> Saved " << out_dir << "/field_drift_velocity_67_33.pdf\n";


        // 4. Mirror Symmetry Diagnostic
        double sum_diff = 0.0, sum_mag = 0.0;
        int nPts = 0;
        for (double yy = -15.9; yy <= 15.9; yy += 0.2) {
            for (double xx = 0.05; xx <= 1.40; xx += 0.05) {   // was 1.49
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
        std::cout << "-> Mirror Symmetry Diagnostic: " << 100.0 * (1.0 - (sum_diff / sum_mag)) << "%\n";
       
        // 5. Electric Field Uniformity Diagnostic
        std::cout << "\n--- FIELD UNIFORMITY DIAGNOSTIC ---\n";
        double sum_e = 0.0;
        double sum_e2 = 0.0;
        int n_uniform_pts = 0;


        for (double x = -1.4; x <= 1.4; x += 0.1) {
            for (double y = -14.0; y <= 14.001; y += 0.2) {
                double ex = 0., ey = 0., ez = 0., v = 0.;
                Medium* m = nullptr;
                int status = 0;
                elm.ElectricField(x, y, 0.005, ex, ey, ez, v, m, status);
               
                if (status == 0) {
                    double emag = std::sqrt(ex*ex + ey*ey + ez*ez);
                    sum_e += emag;
                    sum_e2 += (emag * emag);
                    n_uniform_pts++;
                }
            }
        }


        if (n_uniform_pts > 1) {
            double avg_e = sum_e / n_uniform_pts;
            double variance = (sum_e2 - ((sum_e * sum_e) / n_uniform_pts)) / (n_uniform_pts - 1);
            if (variance < 0.0) variance = 0.0;
            double std_dev = std::sqrt(variance);


            std::cout << "Sampled " << n_uniform_pts << " valid points inside X:[-1.49, 1.49], Y:[-15.0, 15.0]\n";
            std::cout << "Average Field (|E|) : " << avg_e << " V/cm\n";
            std::cout << "Standard Deviation  : " << std_dev << " V/cm\n";
            std::cout << "Field Fluctuation   : " << (std_dev / avg_e) * 100.0 << " %\n";
        }
        std::cout << "-----------------------------------\n\n";


    } else if (switch_mode == 1) {
        std::cout << "\n=========================================\n";
        std::cout << "MODE 1: SPATIAL EFFICIENCY MULTI-GAS SWEEP\n";
        std::cout << "=========================================\n";


        // ==========================================================
        // 1. CALCULATE AVERAGE ELECTRIC FIELD
        // ==========================================================
        double sum_e = 0.0;
        int n_uniform_pts = 0;
        for (double x = -1.4; x <= 1.401; x += 0.1) {   // was -1.49: inside the wire
                                                        // concentration zone (4x at 100um)
            for (double y = -15.0; y <= 15.0; y += 0.2) {
                double ex = 0., ey = 0., ez = 0., v = 0.;
                Medium* m = nullptr;
                int status = 0;
                elm.ElectricField(x, y, 0.005, ex, ey, ez, v, m, status);
               
                if (status == 0) {
                    sum_e += std::sqrt(ex*ex + ey*ey + ez*ez);
                    n_uniform_pts++;
                }
            }
        }
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


        TMultiGraph* mg = new TMultiGraph();
        mg->SetTitle("Total Spatial Efficiency vs Spawn Y;Spawn Y Position [cm];Electrons Passed [%]");
       
        TMultiGraph* mgCurves = new TMultiGraph();
        mgCurves->SetTitle("Spatial Efficiency Fits;Spawn Y Position [cm];Electrons Passed [%]");
       
        auto legend = new TLegend(0.80, 0.52, 0.985, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.035);


        TGraph* grDiff = new TGraph();
        grDiff->SetTitle(Form("Transverse Diffusion at %.1f V/cm;CF4 Concentration [%%];Transverse Diffusion [cm^{1/2}]", avg_e));
        grDiff->SetMarkerStyle(20); grDiff->SetMarkerColor(kBlue); grDiff->SetLineColor(kBlue); grDiff->SetLineWidth(2);


        TGraph* grVel = new TGraph();
        grVel->SetTitle(Form("Drift velocity vs cf4/isobutane concentration at %.1f V/cm;CF4 Concentration [%%];Drift Velocity [cm/ns]", avg_e));
        grVel->SetMarkerStyle(20); grVel->SetMarkerColor(kMagenta+2); grVel->SetLineColor(kMagenta+2); grVel->SetLineWidth(2);


        int graph_idx = 0;


        // ==========================================================
        // 3. RUN GAS SWEEP
        // ==========================================================
        for (int i = 0; i < (int)o2_tags.size(); ++i) {
            const std::string tag   = o2_tags[i];
            const std::string label = o2_labels[i];
            std::string gas_file = "gas_tables_v2/" + tag + ".gas";
            gas.LoadGasFile(gas_file);


            std::string out_name = out_dir + "/epass_" + tag + ".txt";
            std::ofstream outFile(out_name);
            outFile << "Spawn_Y\tElectrons_Passed_Pct\n";
            std::string lost_name = out_dir + "/lost_" + tag + ".txt";
            std::ofstream lostFile(lost_name);
            lostFile << "Spawn_Y\tx0\tx_end\ty_end\tt_end\tnpts\n";
            std::cout << "\nTesting Gas: 67/33 CF4/iC4H10 + O2 " << label << "\n";


            // Record Transverse Diffusion at avg_e
            double dl = 0, dt = 0;
            gas.ElectronDiffusion(avg_e, 0., 0., 0., 0., 0., dl, dt);
            grDiff->SetPoint(graph_idx, i, dt);


            // Record Drift Velocity at avg_e
            double vx = 0., vy = 0., vz = 0.;
            gas.ElectronVelocity(avg_e, 0., 0., 0., 0., 0., vx, vy, vz);
            double v_mag = std::sqrt(vx*vx + vy*vy + vz*vz);
            grVel->SetPoint(graph_idx, i, v_mag);


            TGraph* grSpatial = new TGraph();
            int pt_idx = 0;
            // Derived from y_step so the sweep always spans exactly
            // Y_MIN..Y_MAX inclusive, whatever step size is set above.
            // lround (not truncation) so 0.3 -> 101 points, not 100.
            const int nY = int(std::lround((Y_MAX - Y_MIN) / y_step)) + 1;
            for (int iy = 0; iy < nY; ++iy) {
                const double y_spawn = Y_MIN + y_step * iy;
                const int nElectrons = 500;
                int leaked = 0;
               
                // Serial drift (no OpenMP). Mirrors strip_cage.cc mode 1,
                // including the lost-electron record, so the two cages are parallel.
                for (int e = 0; e < nElectrons; ++e) {
                    const double x0 = -1.4 + (2.8 * e) / (nElectrons - 1.0);
                    drift.DriftElectron(x0, y_spawn, 0.005, 0.0);

                    const int nPoints = drift.GetNumberOfDriftLinePoints();
                    double x1 = 0, y1 = 0, z1 = 0, t1 = 0;
                    if (nPoints > 0) drift.GetDriftLinePoint(nPoints - 1, x1, y1, z1, t1);

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
                std::cout << "  Y = " << y_spawn << " -> " << pass_pct << "% passed\n";
            }


            outFile.close();
            lostFile.close();


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
            std::cout << "-> Saved Individual Plot: " << out_dir << "/" << tag << "_EPass_v_Y.pdf\n";
            delete cSingle;


            TGraph* grCurveOnly = new TGraph(fitFunc);
            grCurveOnly->SetLineColor(colors[graph_idx]);
            grCurveOnly->SetLineWidth(2);
            mgCurves->Add(grCurveOnly, "L");
           
            graph_idx++;
        }


        // ==========================================================
        // 4. FINAL OUTPUTS
        // ==========================================================


        // Spatial Efficiency TOTAL
        TCanvas* cEff = new TCanvas("cEff", "Total Spatial Efficiency", 1000, 600);
        cEff->SetRightMargin(0.22);
        cEff->SetGrid();
        mg->SetMinimum(0.0);
        mg->SetMaximum(100.0);
        mg->Draw("A");
        // GetXaxis() returns nullptr until the multigraph has been drawn, so
        // the X range has to be pinned AFTER Draw("A"), then the pad updated.
        // SetMinimum/SetMaximum above handle Y and work either side of Draw.
        mg->GetXaxis()->SetLimits(Y_MIN, Y_MAX);
        mg->GetXaxis()->SetNdivisions(510);
        gPad->Modified();
        gPad->Update();
        legend->Draw();
        cEff->SaveAs((out_dir + "/Total_EPass_v_Y.pdf").c_str());
        std::cout << "\n-> Saved " << out_dir << "/Total_EPass_v_Y.pdf\n";


        // Spatial Efficiency TOTAL CURVES
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
        std::cout << "-> Saved " << out_dir << "/TotalCurves_EPass_v_Y.pdf\n";


        // Transverse Diffusion
        TCanvas* cDiff = new TCanvas("cDiff", "Transverse Diffusion", 800, 600);
        cDiff->SetGrid();
        grDiff->Draw("APL");
        cDiff->SaveAs((out_dir + "/field_transverse_diffusion_vs_cf4.pdf").c_str());
        std::cout << "-> Saved " << out_dir << "/field_transverse_diffusion_vs_cf4.pdf\n";


        // Drift Velocity
        TCanvas* cVelSweep = new TCanvas("cVelSweep", "Drift Velocity vs CF4", 800, 600);
        cVelSweep->SetGrid();
        grVel->Draw("APL");
        cVelSweep->SaveAs((out_dir + "/field_drift_velocity_vs_cf4.pdf").c_str());
        std::cout << "-> Saved " << out_dir << "/field_drift_velocity_vs_cf4.pdf\n";
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

        std::ofstream fcsv(out_dir + "/efield_full_wire.csv");
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
                  << out_dir << "/efield_full_wire.csv\n";
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
        // One CSV per square: efield_micro_wire_c{C}.csv
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
            fn << out_dir << "/efield_micro_wire_c" << col << ".csv";
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
// =====================================================================
//  field_cage.cc  --  MODE 0 seeding patch
//  Replaces the random-electron block (lines 92-109) with a deterministic
//  n-per-layer grid stepped in y, matching MODE 1's seeding so the two
//  modes are directly comparable and the drift lines are readable.
//  Needs <iomanip> only if you re-add setw(); as written it does not.