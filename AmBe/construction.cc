#include "construction.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
    fMessenger = new G4GenericMessenger(this, "/detector/", "Detector Construction");

    fMessenger->DeclareProperty("nCols", nCols, "Number of columns");
    fMessenger->DeclareProperty("nRows", nRows, "Number of rows");

    fMessenger = new G4GenericMessenger(this, "/germanium/", "Germanium Size");

    fMessenger->DeclareProperty("gerscale", gerscale, "Germanium scale");

    nCols = 100;
    nRows = 100;
    gerscale = 1;
    DefineMaterials();
}

MyDetectorConstruction::~MyDetectorConstruction()
{
}

void MyDetectorConstruction::DefineMaterials()
{
    G4NistManager *nist = G4NistManager::Instance();
    G4double A; //atomic mass
    G4double Z; // atomic number
    G4double d; // density
    A  =  54.94*g/mole;
    G4Element* elMn   =  new G4Element("Manganese","Mn",Z = 25.,A);

    A = 28.09*g/mole;
    G4Element* elSi  = new G4Element("Silicon","Si",Z = 14.,A);

    A = 52.00*g/mole;
    G4Element* elCr  = new G4Element("Chromium","Cr",Z = 24.,A);

    A = 58.70*g/mole;
    G4Element* elNi  = new G4Element("Nickel","Ni",Z = 28.,A);

    A = 55.85*g/mole;
    G4Element* elFe  = new G4Element("Iron","Fe",Z = 26.,A);
  
    worldMat = nist->FindOrBuildMaterial("G4_AIR");

    water = nist->FindOrBuildMaterial("G4_WATER");

    matparaffin = nist->FindOrBuildMaterial("G4_PARAFFIN");


    d = 8.02*g/cm3 ;

    matsteel = new G4Material("Stainless steel",d,5);
    matsteel->AddElement(elMn, 0.02);
    matsteel->AddElement(elSi, 0.01);
    matsteel->AddElement(elCr, 0.19);
    matsteel->AddElement(elNi, 0.10);
    matsteel->AddElement(elFe, 0.68);
    
    /* Creates Materials*/
    fH2O = new G4Material("H2O", 1.000*g/cm3, 2);
    fH2O->AddElement(nist->FindOrBuildElement("H"), 2); 
    fH2O->AddElement(nist->FindOrBuildElement("O"), 1);

    H  = new G4Element("Hydrogen", "H", 1., 1.00794 * g / mole);
    C  = new G4Element("Carbon", "C", 6., 12.011 * g / mole);
    O  = new G4Element("Oxygen", "O", 8., 16.00 * g / mole);
    fWLSfoilPMMA = new G4Material("PMMA", 1.18 * g / cm3, 3);
    fWLSfoilPMMA->AddElement(H, 0.08);
    fWLSfoilPMMA->AddElement(C, 0.60);
    fWLSfoilPMMA->AddElement(O, 0.32);
   
    /* Or use with natural abundances*/
    elGd = nist->FindOrBuildElement("Gd");
    Gdmat = nist->FindOrBuildMaterial("G4_Gd");


    Gdsol = new G4Material("Gd_Solution", 1.000*g/cm3, 2);
    Gdsol->AddMaterial(fH2O, 99.8*perCent);  /*Final Material consisting of above stuff and their %*/
    Gdsol->AddMaterial(Gdmat, 0.2*perCent);



    fworldMat = nist->FindOrBuildMaterial("G4_AIR"); /*Use a Material of the Nistmanager*/

    /*add refractive index of proton in aerogel*/
    fenergySmall = { 1.239841939*eV/0.6,1.239841939*eV/0.25 };
    std::vector<G4double> rindexWorld = {1.0,1.0};
    std::vector<G4double> rindexWater = {1.33,1.33};
    std::vector<G4double> AbslengthWater = {20 * m,20 * m};


    G4MaterialPropertiesTable *mptWater = new G4MaterialPropertiesTable();
    mptWater->AddProperty("RINDEX", fenergySmall, rindexWater, 2); // Thing to add, Momentum, rindex, Anzahl an Wertepaaren
    mptWater->AddProperty("ABSLENGTH", fenergySmall, AbslengthWater, 2);
    fH2O->SetMaterialPropertiesTable(mptWater);

    G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();

    mptWorld->AddProperty("RINDEX", fenergySmall, rindexWorld, 2);

    fworldMat->SetMaterialPropertiesTable(mptWorld);

    fenergy = {
    1.239841939*eV/0.6, 1.239841939*eV/0.59, 1.239841939*eV/0.58, 1.239841939*eV/0.57, 1.239841939*eV/0.56,
    1.239841939*eV/0.55, 1.239841939*eV/0.54, 1.239841939*eV/0.53, 1.239841939*eV/0.52, 1.239841939*eV/0.51,
    1.239841939*eV/0.50, 1.239841939*eV/0.49, 1.239841939*eV/0.48, 1.239841939*eV/0.47, 1.239841939*eV/0.46,
    1.239841939*eV/0.45, 1.239841939*eV/0.44, 1.239841939*eV/0.43, 1.239841939*eV/0.42, 1.239841939*eV/0.41,
    1.239841939*eV/0.40, 1.239841939*eV/0.39, 1.239841939*eV/0.38, 1.239841939*eV/0.37, 1.239841939*eV/0.36,
    1.239841939*eV/0.35, 1.239841939*eV/0.34, 1.239841939*eV/0.33, 1.239841939*eV/0.32, 1.239841939*eV/0.31,
    1.239841939*eV/0.30, 1.239841939*eV/0.29, 1.239841939*eV/0.28, 1.239841939*eV/0.27, 1.239841939*eV/0.26,
    1.239841939*eV/0.25
    };

    std::vector<G4double> refractiveIndexWLSfiber = { 1.60, 1.60 };

    std::vector<G4double> absfiber = {
    0.1 * m, 0.1 *  m, 0.1 * m, 0.1 * m, 0.1 * m, 
    0.1 * m, 0.1 * m, 0.1 * m, 0.1 * m, 0.1 * m,
    0.1 * m, 0.1 * m, 0.1 * m, 0.1 * m, 0.1 * m, 
    0.1 * m, 0.1 * m, 0.1 * m, 0.1 * m, 0.1 * m,
    3 * cm, 0.8 * cm, 2. * mm, 0.09 * mm, 0.08 * mm,
    0.08 * mm, 0.08 * mm,  0.09 * mm,  0.09 * mm,  0.1 * mm, 
    0.12 * mm,  0.09 * mm,  0.08 * mm,  0.08 * mm,  0.11 * mm,
    0.13 * mm
    };

    std::vector<G4double> absWLSfiber = {
    0.08 * m, 0.08 * m, 0.08 * m, 0.08 * m, 0.08 * m, 
    0.08 * m, 0.08 * m, 0.08 * m, 0.08 * m, 0.08 * m,
    0.08 * m, 0.08 * m, 0.08 * m, 0.08 * m, 0.08 * m, 
    0.08 * m, 0.08 * m, 0.08 * m, 0.08 * m, 0.08 * m,
    0.08 * m, 7. * cm, 9.0 * mm, 1.5 * mm, 1.6 * mm,
    1.7 * mm, 1.8 * mm,  1.6 * mm,  1.9 * mm,  2.3 * mm, 
    2.4 * mm,  2.4 * mm,  2.4 * mm,  2.2 * mm,  3. * mm,
    4.4 * mm
    };

    std::vector<G4double> emissionFib = {
    0.05, 0.10, 0.30, 0.50, 0.75, 
    1.00, 1.50, 1.85, 2.30, 2.75,
    3.25, 3.80, 4.50, 5.80, 8.00, 
    10.00, 15.00, 17.0, 12.00, 15.0,
    9.00, 2.50, 1.00, 0.05, 0.00, 
    0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 
    0.00, 
    };

    // Add entries into properties table
    G4MaterialPropertiesTable* mptWLSfiber = new G4MaterialPropertiesTable();
    mptWLSfiber->AddProperty("RINDEX", fenergySmall, refractiveIndexWLSfiber);
    mptWLSfiber->AddProperty("ABSLENGTH", fenergy, absfiber);
    mptWLSfiber->AddProperty("WLSABSLENGTH", fenergy, absWLSfiber);
    mptWLSfiber->AddProperty("WLSCOMPONENT", fenergy, emissionFib);
    mptWLSfiber->AddConstProperty("WLSTIMECONSTANT", 0.5 * ns);

    fWLSfoilPMMA->SetMaterialPropertiesTable(mptWLSfiber);
}

G4Transform3D MyDetectorConstruction::Rotation(G4double theta, G4double x_1, G4double y_1, G4double z_1, G4double x_2, G4double y_2, G4double z_2)
{
    G4Transform3D Trans = G4Rotate3D(theta, G4ThreeVector(x_1, y_1, z_1)) * G4Translate3D(x_2, y_2, z_2);
 
    return Trans;
}
G4Transform3D MyDetectorConstruction::doubleRotation(G4double theta, G4double x_1, G4double y_1, G4double z_1, G4double x_2, G4double y_2, G4double z_2)
{
    G4Transform3D Trans = G4Rotate3D(theta, G4ThreeVector(x_1, y_1, z_1)) * G4Rotate3D(180 * degree, G4ThreeVector(0, 1, 0)) * G4Translate3D(x_2, y_2, z_2);

    return Trans;
}
G4Transform3D MyDetectorConstruction::rotZ(G4double theta, G4double x_1, G4double y_1, G4double z_1)
{
    G4Transform3D Trans = G4RotateZ3D(theta) * G4Translate3D(x_1,y_1,z_1);
 
    return Trans;
}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
    // test//
    G4double alpha = 72 * M_PI / 180;
    G4double beta = 54 * M_PI / 180;
    G4double aWorld = 0.450 * m;                         // mm
    G4double lWorld = aWorld / (2 * sin(36 * degree)); // mm

    G4double r_i = (aWorld / 20) * sqrt(250 + 110 * sqrt(5)); // mm
    G4double r_k = (aWorld / 4) * (3 + sqrt(5));              // mm
    G4double r_e = (aWorld / 4) * sqrt(3) * (1 + sqrt(5));    // mm
    G4double d_l = aWorld / (2 * tan(36 * degree));           // mm

    G4double R_i = (1.113516364) * aWorld; // mm
    G4double R_k = 1.309016994 * aWorld;   // mm

    G4double Beta = acos(r_i / r_k);                                // Rad
    G4double Alpha = acos(r_i / r_e);                               // Rad
    G4double Gamma_div2 = asin(aWorld / (2 * r_e));                 // Rad
    G4double Theta = M_PI - 2 * atan((sqrt(5) - 1) / 2);            // Rad
    G4double Phi = (M_PI / 2) + atan(((sqrt(5) - 1) / 2));          // Rad
    G4double Kappa = 2 * Beta - acos(d_l * cos(36 * degree / r_i)); // Rad

    G4double SteelThickness = 10 * mm ;

    G4double FoilThickness =  65 * um ;//


    //World Cube
    G4double world = 2 * m; // mm
    
    //DoDi Steel
    G4double phiStart = 0;
    G4double phiTotal = 2 * M_PI;
    G4int numSide = 5;
    G4int numZPlanes = 2;
    G4double rInnerSteel[] = {0, 0};
    G4double rOuterSteel[] = {0, (d_l/r_i)*(r_i + SteelThickness + FoilThickness)};
    G4double zPlaneSteel[] = {0, r_i + SteelThickness + FoilThickness};
    
    //DoDi Foil
    G4double rInnerFoil[] = {0, 0};
    G4double rOuterFoil[] = {0, (d_l/r_i)*(r_i+FoilThickness)};
    G4double zPlaneFoil[] = {0, r_i + FoilThickness};

    //DoDi Water
    G4double rInnerWater[] = {0, 0};
    G4double rOuterWater[] = {0, d_l};
    G4double zPlaneWater[] = {0, r_i};
    
    //Gadolinium
    G4double rInnerGadolinium[] = {0, 0};
    G4double rOuterGadolinium[] = {(d_l/r_i)*(r_i - 5*cm), (d_l/r_i)*(r_i -5*cm)};
    G4double zPlaneGadolinium[] = {0,5*cm};

    //PMT Polycone
    G4double realRadius =(252/2) * mm;
    G4double realHeight =-80*mm ;
    G4double numZPlanes1=11;
    G4double zPlane1[]={0,realHeight*0.1,realHeight*0.2,realHeight*0.3,realHeight*0.4,realHeight*0.5,realHeight*0.6,realHeight*0.7,realHeight*0.8,realHeight*0.9,realHeight};
    G4double rInner1[]={0,0,0,0,0,0,0,0,0,0,0};
    G4double rOuter1[]={realRadius,realRadius*0.98,realRadius*0.95,realRadius*0.9,realRadius*0.85,realRadius*0.7,realRadius*0.60,realRadius*0.5,realRadius*0.4,realRadius*0.25,0};
    
    //Paraffin 
    G4double paraffinlength = 5*cm ;

    G4cout << "a: " << aWorld << G4endl;
    G4cout << "l: " << lWorld << G4endl;
    G4cout << "d_l: " << d_l << G4endl;
    G4cout << "r_i/aWorld: " << r_i / aWorld << G4endl;
    G4cout << "r_i: " << r_i << G4endl;
    G4cout << "y: " << r_i * tan(Beta) * sin(72 * degree) << G4endl;
    G4cout << "x: " << -r_i * tan(Beta) * cos(72 * degree) << G4endl;
    G4cout << "r_k/aWorld: " << r_k / aWorld << G4endl;
    G4cout << "beta: " << Beta * 180 / M_PI << G4endl;
    G4cout << "theta: " << Theta * 180 / M_PI << G4endl;
    G4cout << "phi: " << Phi * 180 / M_PI << G4endl;
    G4cout << " pi?:" << 2 * (Alpha + Beta + Gamma_div2) << G4endl;
    G4cout<< "angle: "<< tan(d_l/r_i)*(180/M_PI)<<G4endl;
    G4cout<< "r_k/2m: "<< (r_k/2 * m)<< "tan()"<< atan(r_k/(2 * m)) * (180/ M_PI) << G4endl;


    //solid volume
    solidWorld = new G4Box("solidWorld", world, world, world);
    
    solidSteel = new G4Polyhedra("solidFoil", phiStart, phiTotal, numSide, numZPlanes, zPlaneSteel, rInnerSteel, rOuterSteel);
    solidSteel1 = new G4Polyhedra("solidFoil1", phiStart, phiTotal, numSide, numZPlanes, zPlaneSteel, rInnerSteel, rOuterSteel);

    solidFoil = new G4Polyhedra("solidFoil", phiStart, phiTotal, numSide, numZPlanes, zPlaneFoil, rInnerFoil, rOuterFoil);
    solidFoil1 = new G4Polyhedra("solidFoil1", phiStart, phiTotal, numSide, numZPlanes, zPlaneFoil, rInnerFoil, rOuterFoil);

    solidWater = new G4Polyhedra("solidWater", phiStart, phiTotal, numSide, numZPlanes, zPlaneWater, rInnerWater, rOuterWater);
    solidWater1 = new G4Polyhedra("solidWater1", phiStart, phiTotal, numSide, numZPlanes, zPlaneWater, rInnerWater, rOuterWater);

    solidDetector = new G4Polycone("solidDetector",phiStart, phiTotal, numZPlanes1,zPlane1, rInner1, rOuter1); 
    
    solidPanel = new G4Box("solidPanel", 0.5*m, 0.5*m, 0.01*m);

    solidHuman = new G4Tubs("solidHuman", 0, 40*cm, 1*m, 0, 2*M_PI);

    solidGadolinium = new G4Polyhedra("solidGadolinium", phiStart, phiTotal, numSide, numZPlanes, zPlaneGadolinium, rInnerGadolinium, rOuterGadolinium);

    solidParaffin = new G4Box("solidParaffin", 37.5*cm, 15*cm, 50*cm);
    // logic volume
    logicWorld = new G4LogicalVolume(solidWorld, fworldMat, "logicWorld");

    logicSteel= new G4LogicalVolume(solidSteel, matsteel , "logicSteel");
    logicSteel1= new G4LogicalVolume(solidSteel1, matsteel , "logicSteel");

    logicFoil = new G4LogicalVolume(solidFoil, fWLSfoilPMMA , "logicFoil");
    logicFoil1 = new G4LogicalVolume(solidFoil1, fWLSfoilPMMA , "logicFoil");

    logicWater = new G4LogicalVolume(solidWater, fH2O, "logicWater");
    logicWater1 = new G4LogicalVolume(solidWater1, fH2O, "logicWater");

    logicDetector = new G4LogicalVolume(solidDetector, fworldMat, "logicDetector");

    logicPanel = new G4LogicalVolume(solidPanel, fworldMat, "logicPanel");

    logicHuman = new G4LogicalVolume(solidHuman, fworldMat, "logicHuman");

    logicGadolinium =new G4LogicalVolume(solidGadolinium, Gdsol, "logicGadolinium");

    logicParaffin = new G4LogicalVolume(solidParaffin, matparaffin, "logicParaffin");
   
    //World Vol. placement
    physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);

    //Steel DoDi
    physSteel1 = new G4PVPlacement(
        Rotation(0, 0, 0, r_i+SteelThickness+FoilThickness, 0, 0, 0),
        logicSteel1, "physSteel", logicWorld, false,01, true);
    physSteel2 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+(SteelThickness+FoilThickness)) * tan(Beta), 0, r_i+(SteelThickness+FoilThickness), 0, 0, 0),
        logicSteel, "physSteel", logicWorld, false, 0, true);
    physSteel3 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+(SteelThickness+FoilThickness))* tan(Beta) * cos(72 * degree), (r_i+(SteelThickness+FoilThickness)) * tan(Beta) * sin(72 * degree), r_i+(SteelThickness+FoilThickness), 0, 0, 0),
        logicSteel, "physSteel", logicWorld, false, 0, true);
    physSteel4 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+(SteelThickness+FoilThickness)) * tan(Beta) * cos(2 * 72 * degree), (r_i+(SteelThickness+FoilThickness)) * tan(Beta) * sin(2 * 72 * degree), r_i+(SteelThickness+FoilThickness), 0, 0, 0),
        logicSteel, "physSteel", logicWorld, false, 0, true);
    physSteel5 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+(SteelThickness+FoilThickness)) * tan(Beta) * cos(72 * degree), -(r_i+(SteelThickness+FoilThickness)) * tan(Beta) * sin(72 * degree), r_i+(SteelThickness+FoilThickness), 0, 0, 0),
        logicSteel, "physSteel", logicWorld, false, 0, true);
    physSteel6 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+(SteelThickness+FoilThickness)) * tan(Beta) * cos(2 * 72 * degree), -(r_i+(SteelThickness+FoilThickness)) * tan(Beta) * sin(2 * 72 * degree), r_i+(SteelThickness+FoilThickness), 0, 0, 0),
        logicSteel, "physSteel", logicWorld, false, 0, true);
    physSteel7 = new G4PVPlacement(
        doubleRotation(0, 0, 0, r_i+(SteelThickness+FoilThickness), 0, 0, 0),
        logicSteel, "physSteel", logicWorld, false, 0, true);
    physSteel8 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+(SteelThickness+FoilThickness)) * tan(Beta), 0, r_i+(SteelThickness+FoilThickness), 0, 0, 0),
        logicSteel, "physSteel", logicWorld, false, 0, true);
    physSteel9 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+(SteelThickness+FoilThickness)) * tan(Beta) * cos(72 * degree), (r_i+(SteelThickness+FoilThickness)) * tan(Beta) * sin(72 * degree), r_i+(SteelThickness+FoilThickness), 0, 0, 0),
        logicSteel, "physSteel", logicWorld, false, 0, true);
    physSteel10 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+(SteelThickness+FoilThickness)) * tan(Beta) * cos(2 * 72 * degree), (r_i+(SteelThickness+FoilThickness)) * tan(Beta) * sin(2 * 72 * degree), r_i+(SteelThickness+FoilThickness), 0, 0, 0),
        logicSteel, "physSteel", logicWorld, false, 0, true);
    physSteel11 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+(SteelThickness+FoilThickness)) * tan(Beta) * cos(72 * degree), -(r_i+(SteelThickness+FoilThickness)) * tan(Beta) * sin(72 * degree), r_i+(SteelThickness+FoilThickness), 0, 0, 0),
        logicSteel, "physSteel", logicWorld, false, 0, true);
    physSteel12 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+(SteelThickness+FoilThickness)) * tan(Beta) * cos(2 * 72 * degree), -(r_i+(SteelThickness+FoilThickness)) * tan(Beta) * sin(2 * 72 * degree), r_i+(SteelThickness+FoilThickness), 0, 0, 0),
        logicSteel, "physSteel", logicWorld, false, 0, true);

    //Foil DoDi
    physFoil1 = new G4PVPlacement(
        Rotation(0, 0, 0, r_i+FoilThickness, 0, 0, 0),
        logicFoil, "physFoil", logicSteel,  false,0, true);
    /*physFoil2 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+FoilThickness) * tan(Beta), 0, r_i+FoilThickness, 0, 0, 0),
        logicFoil, "physFoil", logicWorld, false, 0, true);
    physFoil3 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+FoilThickness)* tan(Beta) * cos(72 * degree), (r_i+FoilThickness) * tan(Beta) * sin(72 * degree), r_i+FoilThickness, 0, 0, 0),
        logicFoil, "physFoil", logicWorld, false, 0, true);
    physFoil4 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+FoilThickness) * tan(Beta) * cos(2 * 72 * degree), (r_i+FoilThickness) * tan(Beta) * sin(2 * 72 * degree), r_i+FoilThickness, 0, 0, 0),
        logicFoil, "physFoil", logicWorld, false, 0, true);
    physFoil5 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+FoilThickness) * tan(Beta) * cos(72 * degree), -(r_i+FoilThickness) * tan(Beta) * sin(72 * degree), r_i+FoilThickness, 0, 0, 0),
        logicFoil, "physFoil", logicWorld, false, 0, true);
    physFoil6 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+FoilThickness) * tan(Beta) * cos(2 * 72 * degree), -(r_i+FoilThickness) * tan(Beta) * sin(2 * 72 * degree), r_i+FoilThickness, 0, 0, 0),
        logicFoil, "physFoil", logicWorld, false, 0, true);*/
    physFoil7 = new G4PVPlacement(
        Rotation(0, 0, 0, r_i+FoilThickness, 0, 0, 0),
        logicFoil1, "physFoil", logicSteel1, false, 01, true);
/*    physFoil8 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+FoilThickness) * tan(Beta), 0, r_i+FoilThickness, 0, 0, 0),
        logicFoil, "physFoil", logicWorld, false, 0, true);
    physFoil9 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+FoilThickness) * tan(Beta) * cos(72 * degree), (r_i+FoilThickness) * tan(Beta) * sin(72 * degree), r_i+FoilThickness, 0, 0, 0),
        logicFoil, "physFoil", logicWorld, false, 0, true);
    physFoil10 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+FoilThickness) * tan(Beta) * cos(2 * 72 * degree), (r_i+FoilThickness) * tan(Beta) * sin(2 * 72 * degree), r_i+FoilThickness, 0, 0, 0),
        logicFoil, "physFoil", logicWorld, false, 0, true);
    physFoil11 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+FoilThickness) * tan(Beta) * cos(72 * degree), -(r_i+FoilThickness) * tan(Beta) * sin(72 * degree), r_i+FoilThickness, 0, 0, 0),
        logicFoil, "physFoil", logicWorld, false, 0, true);
    physFoil12 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+FoilThickness) * tan(Beta) * cos(2 * 72 * degree), -(r_i+FoilThickness) * tan(Beta) * sin(2 * 72 * degree), r_i+FoilThickness, 0, 0, 0),
        logicFoil, "physFoil", logicWorld, false, 0, true);*/

    //DoDi Water
    physWater1 = new G4PVPlacement(
        Rotation(0, 0, 0, r_i, 0, 0, 0),
        logicWater, "physWater", logicFoil, false, 0, true);
    /*physWater2 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i) * tan(Beta), 0, r_i, 0, 0, 0),
        logicWater, "physWater", logicWorld, false, 0, true);
    physWater3 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i)* tan(Beta) * cos(72 * degree), (r_i) * tan(Beta) * sin(72 * degree), r_i, 0, 0, 0),
        logicWater, "physWater", logicWorld, false, 0, true);
    physWater4 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i) * tan(Beta) * cos(2 * 72 * degree), (r_i) * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, 0),
        logicWater, "physWater", logicWorld, false, 0, true);
    physWater5 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i) * tan(Beta) * cos(72 * degree), -(r_i) * tan(Beta) * sin(72 * degree), r_i, 0, 0, 0),
        logicWater, "physWater", logicWorld, false, 0, true);
    physWater6 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i) * tan(Beta) * cos(2 * 72 * degree), -(r_i) * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, 0),
        logicWater, "physWater", logicWorld, false, 0, true);*/
    physWater7 = new G4PVPlacement(
        Rotation(0, 0, 0, r_i, 0, 0, 0),
        logicWater1, "physWater", logicFoil1, false, 01, true);
/*    physWater8 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i) * tan(Beta), 0, r_i, 0, 0, 0),
        logicWater, "physWater", logicWorld, false, 0, true);
    physWater9 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i) * tan(Beta) * cos(72 * degree), (r_i) * tan(Beta) * sin(72 * degree), r_i, 0, 0, 0),
        logicWater, "physWater", logicWorld, false, 0, true);
    physWater10 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i) * tan(Beta) * cos(2 * 72 * degree), (r_i) * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, 0),
        logicWater, "physWater", logicWorld, false, 0, true);
    physWater11 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i) * tan(Beta) * cos(72 * degree), -(r_i) * tan(Beta) * sin(72 * degree), r_i, 0, 0, 0),
        logicWater, "physWater", logicWorld, false, 0, true);
    physWater12 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i) * tan(Beta) * cos(2 * 72 * degree), -(r_i) * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, 0),
        logicWater, "physWater", logicWorld, false, 0, true);*/

    //DoDi Detector
    physDetector = new G4PVPlacement(
        Rotation(180 * degree,0, 0, r_i, 0, 0, r_i),
        logicDetector, "physDetector", logicWater, false, 2, true);

    //Paraffin Placement
    //physParaffin1 = new G4PVPlacement(rotZ(0, 0, 22.5*cm,110*cm),logicParaffin,"physParaffin",logicWorld,false,0,true);
    //physParaffin2 = new G4PVPlacement(rotZ(0, 0, -22.5*cm,110*cm),logicParaffin,"physParaffin",logicWorld,false,0,true);
    //physParaffin3 = new G4PVPlacement(rotZ(90*degree, 0, 22.5*cm,110*cm),logicParaffin,"physParaffin",logicWorld,false,0,true);
    //physParaffin4 = new G4PVPlacement(rotZ(90*degree, 0, -22.5*cm,110*cm), logicParaffin,"physParaffin",logicWorld,false,0,true);
    
    //physHuman = new G4PVPlacement(0,G4ThreeVector(0,1.5*m,0),logicHuman, "physHuman", logicWorld, false, 69, true);

    physGadolinium = new G4PVPlacement(0, G4ThreeVector(0,0,+(r_i - 5*cm-2*FoilThickness)), logicGadolinium, "physGadolinium", logicWater1, false, 35, true);

    fScoringVolume = logicGadolinium;
    return physWorld;
}

void MyDetectorConstruction::defineBoundaries()
{   
    //Invisible boundary
    G4MaterialPropertiesTable *IMPT = new G4MaterialPropertiesTable();
    std::vector<G4double> reflectivityI = {0.,0.};
    std::vector<G4double> transmissionI = {1.,1.};
    IMPT->AddProperty("REFLECTIVITY", fenergySmall, reflectivityI, 2);
    IMPT->AddProperty("TRANSMITTANCE", fenergySmall, transmissionI, 2);

    G4OpticalSurface *InvisibleBoundary = new G4OpticalSurface("Invisible");
    InvisibleBoundary->SetMaterialPropertiesTable(IMPT);

    //First boundary
    G4MaterialPropertiesTable *firstMPT = new G4MaterialPropertiesTable();
    std::vector<G4double> reflectivityFirst = {0.09,0.09};
    std::vector<G4double> transmissionFirst = {0.91,0.91};
    firstMPT->AddProperty("REFLECTIVITY", fenergySmall, reflectivityFirst, 2);
    firstMPT->AddProperty("TRANSMITTANCE", fenergySmall, transmissionFirst, 2);

    G4OpticalSurface *FirstBoundary = new G4OpticalSurface("FirstBoundary");
    FirstBoundary->SetType(dielectric_dielectric);
    FirstBoundary->SetModel(unified);
    FirstBoundary->SetFinish(groundfrontpainted);
    FirstBoundary->SetMaterialPropertiesTable(firstMPT);

    //Second boundary
    G4MaterialPropertiesTable *secondMPT = new G4MaterialPropertiesTable();
    std::vector<G4double> reflectivitySecond = {1.0,1.0};
    std::vector<G4double> transmissionSecond = {0.,0.};
    secondMPT->AddProperty("REFLECTIVITY", fenergySmall, reflectivitySecond, 2);
    secondMPT->AddProperty("TRANSMITTANCE", fenergySmall, transmissionSecond, 2);

    G4OpticalSurface *SecondBoundary = new G4OpticalSurface("SecondBoundary");
    SecondBoundary->SetType(dielectric_dielectric);
    SecondBoundary->SetModel(unified);
    SecondBoundary->SetFinish(groundfrontpainted);
    SecondBoundary->SetMaterialPropertiesTable(secondMPT);

    G4LogicalBorderSurface *SecondSurface1 = new G4LogicalBorderSurface("SecondSurface", physFoil1, physSteel1 , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface2 = new G4LogicalBorderSurface("SecondSurface", physFoil7, physSteel7 , SecondBoundary);
    /*
    G4LogicalBorderSurface *SecondSurface3 = new G4LogicalBorderSurface("SecondSurface", physFoil3, physSteel3 , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface4 = new G4LogicalBorderSurface("SecondSurface", physFoil4, physSteel4 , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface5 = new G4LogicalBorderSurface("SecondSurface", physFoil5, physSteel5 , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface6 = new G4LogicalBorderSurface("SecondSurface", physFoil6, physSteel6 , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface7 = new G4LogicalBorderSurface("SecondSurface", physFoil7, physSteel7 , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface8 = new G4LogicalBorderSurface("SecondSurface", physFoil8, physSteel8 , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface9 = new G4LogicalBorderSurface("SecondSurface", physFoil9, physSteel9 , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface10 = new G4LogicalBorderSurface("SecondSurface", physFoil10, physSteel10 , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface11 = new G4LogicalBorderSurface("SecondSurface", physFoil11, physSteel11 , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface12 = new G4LogicalBorderSurface("SecondSurface", physFoil12, physSteel12 , SecondBoundary);
    */
}

void MyDetectorConstruction::ConstructSDandField()
{
    MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");

    logicDetector->SetSensitiveDetector(sensDet);
    //logicHuman->SetSensitiveDetector(sensDet);

}
