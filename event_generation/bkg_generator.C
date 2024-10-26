
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <bitset>
#include <TF1.h>
#include <TRandom.h>
#include <filesystem>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <map>

namespace fs = std::filesystem;

// Function to split a string by spaces and return a vector of values
std::vector<std::string> splitLine(const std::string &line)
{
  std::vector<std::string> result;
  std::istringstream iss(line);
  std::string value;
  while (iss >> value)
  {
    result.push_back(value);
  }
  return result;
}

void listFilesInDirectory(const std::string &path)
{
  try
  {
    // Iterate over the directory
    for (const auto &entry : fs::directory_iterator(path))
    {
      // Check if it's a file and print its path
      if (fs::is_regular_file(entry.status()))
      {
        std::cout << entry.path().string() << std::endl;
      }
    }
  }
  catch (const fs::filesystem_error &e)
  {
    std::cerr << "Error accessing directory: " << e.what() << std::endl;
  }
}

float getZ(std::string val){
  float z = 0;
  if(val=="PixStn0")
    z = 71.1750031;
  else if(val=="PixStn1")
    z = 151.175003;
  else if(val=="PixStn2")
    z = 201.175003;
  else if(val=="PixStn3")
    z = 251.175003;
  else if(val=="PixStn4")
    z = 381.174988;
  else if(val=="MS0")
    z = 3000;
  else if(val=="MS1")
    z = 3600;
  else if(val=="MS2")
    z = 5300;
  else if(val=="MS3")
    z = 5900;
  else if(val=="MS4")
    z = 8100;
  else if(val=="MS5")
    z = 8500;
  return z;
}
int getSurface(std::string val){
  int surface = 0;
  if(val=="PixStn0")
    surface = 2;
  else if(val=="PixStn1")
    surface = 4;
  else if(val=="PixStn2")
    surface = 6;
  else if(val=="PixStn3")
    surface = 8;
  else if(val=="PixStn4")
    surface = 10;
  else if(val=="MS0")
    surface = 12;
  else if(val=="MS1")
    surface = 14;
  else if(val=="MS2")
    surface = 16;
  else if(val=="MS3")
    surface = 18;
  else if(val=="MS4")
    surface = 20;
  else if(val=="MS5")
    surface = 22;
  return surface;
}

std::map<int, int> flukaToPdg = {
    // Photon, Electron, Positron
    {7, 22},  // Photon
    {3, 11},  // Electron
    {4, -11}, // Positron

    // Neutrons, Protons, Deuterons, Alphas
    {8, 2112},        // Neutron
    {1, 2212},        // Proton
    {45, 1000010020}, // Deuteron
    {47, 1000020040}, // Alpha particle

    // Charged Pions, Neutral Pions
    {9, 211},   // Charged Pion (π⁺)
    {10, -211}, // Charged Pion (π⁻)
    {13, 111},  // Neutral Pion (π⁰)

    // Muons
    {5, 13},  // Muon (μ⁻)
    {6, -13}, // Muon (μ⁺)

    // Kaons
    {11, 321},  // Kaon (K⁺)
    {12, -321}, // Kaon (K⁻)
    {16, 310},  // Kaon (K⁰_S)
    {15, 130},  // Kaon (K⁰_L)

    // Other hadrons
    {14, 221},  // Eta (η)
    {18, 3122}, // Lambda (Λ⁰)
    {19, 3222}, // Sigma (Σ⁺)
    {20, 3212}, // Sigma (Σ⁰)
    {21, 3112}, // Sigma (Σ⁻)
    {22, 3322}, // Xi (Ξ⁰)
    {23, 3312}, // Xi (Ξ⁻)
    {24, 3334}, // Omega (Ω⁻)

    // Hyperons
    {25, 4112}, // Sigma_c⁺
    {26, 4114}, // Sigma_c*⁺
    {27, 4212}, // Xi_c⁰

    // Baryons and others
    {48, 3124}, // Delta⁰
    {49, 3224}, // Delta⁺

    // Leptons
    {31, 12},  // Neutrino (ν_e)
    {32, -12}, // Anti-Neutrino (ν̅_e)
    {33, 14},  // Muon Neutrino (ν_μ)
    {34, -14}, // Anti-Muon Neutrino (ν̅_μ)

    // Strange particles
    {35, 130},  // Kaon_long (K⁰_L)
    {36, 310},  // Kaon_short (K⁰_S)
    {37, 311},  // Kaon_zero (K⁰)
    {38, -311}, // Anti-Kaon_zero (K̅⁰)

    // Heavy mesons
    {50, 431},  // D_s+
    {51, -431}, // D_s-

    // Antiparticles
    {2, -2212},  // Antiproton
    {17, -3122}, // Anti-Lambda
    {28, -3222}, // Anti-Sigma⁺
    {29, -3112}, // Anti-Sigma⁻
    {30, -3322}, // Anti-Xi⁰
    {46, -3334}, // Anti-Omega⁻
};

void bkg_generator(int nev = -1, int Eint = 40)
{
  /*

  // clang-format off
  /// (2^8)-1 = 255 volumes
  static constexpr Value kVolumeMask    = 0xff00000000000000;
  /// (2^8)-1 = 255 boundaries
  static constexpr Value kBoundaryMask  = 0x00ff000000000000;
  /// (2^12)-1 = 4095 layers
  static constexpr Value kLayerMask     = 0x0000fff000000000;
  /// (2^8)-1 = 255 approach surfaces
  static constexpr Value kApproachMask  = 0x0000000ff0000000;
  static constexpr Value kPassiveMask   = kApproachMask;
  /// (2^20)-1 = 1048575 sensitive surfaces
  static constexpr Value kSensitiveMask = 0x000000000fffff00;
  /// (2^8)-1 = 255 extra values
  static constexpr Value kExtraMask     = 0x00000000000000ff;

  volumes boundaries layers approach-surfaces  sensitive-surfaces extra
  00000001 00000000 000000000100 00000000 00000000000000000001 00000000
        1         0        4        0                        1       0

  */
  float zMS = 2000; // mm
  char csvname_vt[80];
  char csvname_ms[80];
  char csvname_vt_ms[80];

  std::string directoryPath = "inputData/fluka_" + std::to_string(Eint)+ "GeV";

  std::string directoryName_vt = "simulatedEvents/bkghits_" + std::to_string(Eint) + "GeV_vt";
  std::string directoryName_ms = "simulatedEvents/bkghits_" + std::to_string(Eint) + "GeV_ms";
  std::string directoryName_vt_ms = "simulatedEvents/bkghits_" + std::to_string(Eint) + "GeV_vt_ms";

  if (!fs::exists(directoryName_vt))
    bool created = fs::create_directory(directoryName_vt);
    
  if (!fs::exists(directoryName_ms))
    bool created = fs::create_directory(directoryName_ms);
    
  if (!fs::exists(directoryName_vt_ms))
    bool created = fs::create_directory(directoryName_vt_ms);
    
  std::string qa_file_name = directoryName_vt_ms + "/QA_plots.root";
  TFile *check = new TFile(qa_file_name.c_str(), "recreate");

  FILE *fpcsv_vt;
  FILE *fpcsv_ms;
  FILE *fpcsv_vt_ms;

  std::string eventFile = directoryName_vt + "/event" + std::string(8, '0') + std::to_string(static_cast<int>(0)) + "-hits.csv";
  sprintf(csvname_vt, eventFile.c_str(), 0);
  eventFile = directoryName_ms + "/event" + std::string(8, '0') + std::to_string(static_cast<int>(0)) + "-hits.csv";
  sprintf(csvname_ms, eventFile.c_str(), 0);
  eventFile = directoryName_vt_ms + "/event" + std::string(8, '0') + std::to_string(static_cast<int>(0)) + "-hits.csv";
  sprintf(csvname_vt_ms, eventFile.c_str(), 0);

  fpcsv_vt = fopen(csvname_vt, "w");
  fpcsv_ms = fopen(csvname_ms, "w");
  fpcsv_vt_ms = fopen(csvname_vt_ms, "w");

  fprintf(fpcsv_vt, "particle_id,geometry_id,tx,ty,tz,tt,tpx,tpy,tpz,te,deltapx,deltapy,deltapz,deltae,index\n");
  fprintf(fpcsv_ms, "particle_id,geometry_id,tx,ty,tz,tt,tpx,tpy,tpz,te,deltapx,deltapy,deltapz,deltae,index\n");
  fprintf(fpcsv_vt_ms, "particle_id,geometry_id,tx,ty,tz,tt,tpx,tpy,tpz,te,deltapx,deltapy,deltapz,deltae,index\n");
  int bkg_counter = 10;

  int eventNumber = 0;
  try
  {
    // Iterate over the directory
    for (const auto &entry : fs::directory_iterator(directoryPath))
    {
      std::ifstream fileFluka(entry.path().string()); // Replace with your filename
      if (!fileFluka.is_open())
      {
        std::cerr << "Error: Could not open file!" << entry.path().string() << std::endl;
        continue;
      }

      std::cout << entry.path().string() << std::endl;

      std::string line;
      while (std::getline(fileFluka, line))
      {
        // Check if it's the "end of event" line
        if (line.find("**************** end of event **************") != std::string::npos)
        {
          fclose(fpcsv_vt);
          fclose(fpcsv_ms);
          fclose(fpcsv_vt_ms);
          eventNumber++;
          int evTmp = eventNumber;
          if (eventNumber == nev)
          {
            check->Close();
            return;
          }
          int nzeros = 9;
          while (evTmp != 0)
          {
            evTmp /= 10;
            nzeros--;
          }

          std::string eventFile = directoryName_vt + "/event" + std::string(nzeros, '0') + std::to_string(static_cast<int>(eventNumber)) + "-hits.csv";
          sprintf(csvname_vt, eventFile.c_str(), eventNumber);
          eventFile = directoryName_ms + "/event" + std::string(nzeros, '0') + std::to_string(static_cast<int>(eventNumber)) + "-hits.csv";
          sprintf(csvname_ms, eventFile.c_str(), eventNumber);
          eventFile = directoryName_vt_ms + "/event" + std::string(nzeros, '0') + std::to_string(static_cast<int>(eventNumber)) + "-hits.csv";
          sprintf(csvname_vt_ms, eventFile.c_str(), eventNumber);

          fpcsv_vt = fopen(csvname_vt, "w");
          fpcsv_ms = fopen(csvname_ms, "w");
          fpcsv_vt_ms = fopen(csvname_vt_ms, "w");

          fprintf(fpcsv_vt, "particle_id,geometry_id,tx,ty,tz,tt,tpx,tpy,tpz,te,deltapx,deltapy,deltapz,deltae,index\n");
          fprintf(fpcsv_ms, "particle_id,geometry_id,tx,ty,tz,tt,tpx,tpy,tpz,te,deltapx,deltapy,deltapz,deltae,index\n");
          fprintf(fpcsv_vt_ms, "particle_id,geometry_id,tx,ty,tz,tt,tpx,tpy,tpz,te,deltapx,deltapy,deltapz,deltae,index\n");
        }
        else
        {
          // FLUKA format
          // setup element - particle id - int.origin - etot - x - y - z - director cosines
          // PixStn2            3 AbsoPlu1   9.07811627E-04   14.2597504 -10.3486423       20.1179695       3.07518914E-02 6.16198704E-02 0.997625828
          // ACTS hit format
          // particle_id,geometry_id,tx,ty,tz,tt,tpx,tpy,tpz,te,deltapx,deltapy,deltapz,deltae,index
          // 5156621573355995136,72057731476881664,0.307201,-1.999749,3000.000000,0.000000,-182.768341,-211.306107,91.450699,293.969513,0.000000,0.000000,0.000000,0.000000,0
          // tx,ty,tz,tt = hit coordinates
          // tpx,tpy,tpz,te = particle momentum and energy 
          // deltapx,deltapy,deltapz,deltae = variation of momentum and energy
          //
          // always a new particles id

          std::string value;
          std::vector<std::string> lineValues;

          std::istringstream iss(line);
          while (iss >> value)
            lineValues.push_back(value);

          int nValues = lineValues.size();
          int offSet = 0;
          if(nValues==9)
            offSet = -1;
          // Store the parsed values if not empty
          if (!lineValues.empty())
          {
            std::string binaryStringPart = std::bitset<12>(bkg_counter++).to_string() + std::bitset<12>(0).to_string() + std::bitset<16>(1).to_string() + std::bitset<8>(0).to_string() + std::bitset<16>(0).to_string();
            unsigned long long particle_id = std::stoull(binaryStringPart, nullptr, 2);

            //float value = std::stof(str);
            float x = std::stof(lineValues[3+offSet]);
            float y = std::stof(lineValues[4+offSet]);
            //float z = std::stof(lineValues[5+offSet]);
            float e = std::stof(lineValues[6+offSet]);
            float cosx = std::stof(lineValues[7+offSet]);
            float cosy = std::stof(lineValues[8+offSet]);
            float cosz = std::stof(lineValues[9+offSet]);
            int flukacode = std::stof(lineValues[1]);
            float z = getZ(lineValues[0]);
            
            double mass = TDatabasePDG::Instance()->GetParticle(flukaToPdg[flukacode])->Mass();
            
            float p = TMath::Sqrt(e*e-mass*mass);
            float px = p*cosx;
            float py = p*cosy;
            float pz = p*cosz;
            
            int surface = getSurface(lineValues[0]);
            if(surface==0){
              continue;
            }
            int surface_ms = getSurface(lineValues[0])-10;

            if (z < zMS)
            {
              if(flukacode!=3 && flukacode!=4 && z<400){
                continue;
              }
              // volumes (8 bit) + boundaries (8 bit) layers (12 bit) + approach-surfaces (8 bit) + sensitive-surfaces (20 bit) + extra (8 bit)
              std::string binaryStringGeo = std::bitset<8>(1).to_string() + std::string(8, '0') + std::bitset<12>(surface).to_string() + std::string(8, '0') + std::bitset<20>(1).to_string() + std::string(8, '0');
              unsigned long long geometry_id = std::stoull(binaryStringGeo, nullptr, 2);
              fprintf(fpcsv_vt, "%llu,%llu,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%i\n", particle_id, geometry_id, x, y, z, 0., px, py, pz, e, 0., 0., 0., 0., 0);
            }
            else
            {
              std::string binaryStringGeo = std::bitset<8>(1).to_string() + std::string(8, '0') + std::bitset<12>(surface_ms).to_string() + std::string(8, '0') + std::bitset<20>(1).to_string() + std::string(8, '0');
              unsigned long long geometry_id = std::stoull(binaryStringGeo, nullptr, 2);
              fprintf(fpcsv_ms, "%llu,%llu,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%i\n", particle_id, geometry_id, x, y, z, 0., px, py, pz, e, 0., 0., 0., 0., 0);
            }
            std::string binaryStringGeo = std::bitset<8>(1).to_string() + std::string(8, '0') + std::bitset<12>(surface).to_string() + std::string(8, '0') + std::bitset<20>(1).to_string() + std::string(8, '0');
            unsigned long long geometry_id = std::stoull(binaryStringGeo, nullptr, 2);
            fprintf(fpcsv_vt_ms, "%llu,%llu,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%i\n", particle_id, geometry_id, x, y, z, 0., px, py, pz, e, 0., 0., 0., 0., 0);
          }
        }
      }
      fileFluka.close();
    }
  }
  catch (const fs::filesystem_error &e)
  {
    std::cerr << "Error accessing directory: " << e.what() << std::endl;
  }
  check->Close();
}
