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

    worldMat = nist->FindOrBuildMaterial("G4_AIR");
    water = nist->FindOrBuildMaterial("G4_WATER");
    
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

    G4double thickness =  65 * um ;//


    //World Cube
    G4double world = 2 * m; // mm
    
    //DoDi Surf
    G4double phiStart = 0;
    G4double phiTotal = 2 * M_PI;
    G4int numSide = 5;
    G4int numZPlanes = 2;
    G4double rInner0[] = {0, 0};
    G4double rOuter0[] = {0, (d_l/r_i)*(r_i+thickness)};
    G4double zPlane0[] = {0, r_i + thickness};

    //DoDi
    G4double rInner[] = {0, 0};
    G4double rOuter[] = {0, d_l};
    G4double zPlane[] = {0, r_i};

    //PMT Polycone
    G4double realRadius =(252/2) * mm;
    G4double realHeight =-80*mm ;
    G4double numZPlanes1=11;
    G4double zPlane1[]={0,realHeight*0.1,realHeight*0.2,realHeight*0.3,realHeight*0.4,realHeight*0.5,realHeight*0.6,realHeight*0.7,realHeight*0.8,realHeight*0.9,realHeight};
    G4double rInner1[]={0,0,0,0,0,0,0,0,0,0,0};
    G4double rOuter1[]={realRadius,realRadius*0.98,realRadius*0.95,realRadius*0.9,realRadius*0.85,realRadius*0.7,realRadius*0.60,realRadius*0.5,realRadius*0.4,realRadius*0.25,0};


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

    //solid volume
    solidWorld = new G4Box("solidWorld", world, world, world);

    solidSurf = new G4Polyhedra("solidSurf", phiStart, phiTotal, numSide, numZPlanes, zPlane0, rInner0, rOuter0);

    solidDoDi = new G4Polyhedra("solidDoDi", phiStart, phiTotal, numSide, numZPlanes, zPlane, rInner, rOuter);

    solidDetector = new G4Polycone("solidDetector",phiStart, phiTotal, numZPlanes1,zPlane1, rInner1, rOuter1); 
    
    solidPanel = new G4Box("solidPanel", 0.5*m, 0.5*m, 0.01*m);

    // logic volume
    logicWorld = new G4LogicalVolume(solidWorld, fworldMat, "logicWorld");

    logicSurf = new G4LogicalVolume(solidSurf,fWLSfoilPMMA , "logicSurf");

    logicDoDi = new G4LogicalVolume(solidDoDi, fH2O, "logicDoDi");

    logicDetector = new G4LogicalVolume(solidDetector, fworldMat, "logicDetector");

    logicPanel = new G4LogicalVolume(solidPanel, fworldMat, "logicPanel");
   
    //World Vol. placement
    physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);

    //SURF DoDi
    physSurf1 = new G4PVPlacement(
        Rotation(0, 0, 0, r_i+thickness, 0, 0, 0),
        logicSurf, "physSurf", logicWorld, false,111, true);
    physSurf2 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+thickness) * tan(Beta), 0, r_i+thickness, 0, 0, 0),
        logicSurf, "physSurf", logicWorld, false, 112, true);
    physSurf3 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+thickness)* tan(Beta) * cos(72 * degree), (r_i+thickness) * tan(Beta) * sin(72 * degree), r_i+thickness, 0, 0, 0),
        logicSurf, "physSurf", logicWorld, false, 113, true);
    physSurf4 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+thickness) * tan(Beta) * cos(2 * 72 * degree), (r_i+thickness) * tan(Beta) * sin(2 * 72 * degree), r_i+thickness, 0, 0, 0),
        logicSurf, "physSurf", logicWorld, false, 114, true);
    physSurf5 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+thickness) * tan(Beta) * cos(72 * degree), -(r_i+thickness) * tan(Beta) * sin(72 * degree), r_i+thickness, 0, 0, 0),
        logicSurf, "physSurf", logicWorld, false, 115, true);
    physSurf6 = new G4PVPlacement(
        Rotation(180 * degree, -(r_i+thickness) * tan(Beta) * cos(2 * 72 * degree), -(r_i+thickness) * tan(Beta) * sin(2 * 72 * degree), r_i+thickness, 0, 0, 0),
        logicSurf, "physSurf", logicWorld, false, 116, true);
    physSurf7 = new G4PVPlacement(
        doubleRotation(0, 0, 0, r_i+thickness, 0, 0, 0),
        logicSurf, "physSurf", logicWorld, false, 117, true);
    physSurf8 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+thickness) * tan(Beta), 0, r_i+thickness, 0, 0, 0),
        logicSurf, "physSurf", logicWorld, false, 118, true);
    physSurf9 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+thickness) * tan(Beta) * cos(72 * degree), (r_i+thickness) * tan(Beta) * sin(72 * degree), r_i+thickness, 0, 0, 0),
        logicSurf, "physSurf", logicWorld, false, 119, true);
    physSurf10 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+thickness) * tan(Beta) * cos(2 * 72 * degree), (r_i+thickness) * tan(Beta) * sin(2 * 72 * degree), r_i+thickness, 0, 0, 0),
        logicSurf, "physSurf", logicWorld, false, 120, true);
    physSurf11 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+thickness) * tan(Beta) * cos(72 * degree), -(r_i+thickness) * tan(Beta) * sin(72 * degree), r_i+thickness, 0, 0, 0),
        logicSurf, "physSurf", logicWorld, false, 121, true);
    physSurf12 = new G4PVPlacement(
        doubleRotation(180 * degree, -(r_i+thickness) * tan(Beta) * cos(2 * 72 * degree), -(r_i+thickness) * tan(Beta) * sin(2 * 72 * degree), r_i+thickness, 0, 0, 0),
        logicSurf, "physSurf", logicWorld, false, 122, true);
    //DoDi Placement
    physDoDi1 = new G4PVPlacement(
        Rotation(0, 0, 0, r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 1, true);
    physDoDi2 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta), 0, r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 2, true);
  physDoDi3 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(72 * degree), r_i * tan(Beta) * sin(72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 3, true);
    physDoDi4 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 4, true);
    physDoDi5 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(72 * degree), -r_i * tan(Beta) * sin(72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 5, true);
    physDoDi6 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), -r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 6, true);
    physDoDi7 = new G4PVPlacement(
        doubleRotation(0, 0, 0, r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 7, true);
    physDoDi8 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta), 0, r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 8, true);
    physDoDi9 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(72 * degree), r_i * tan(Beta) * sin(72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 9, true);
    physDoDi10 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 10, true);
    physDoDi11 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(72 * degree), -r_i * tan(Beta) * sin(72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 11, true);
    physDoDi12 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), -r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 12, true);

    //Detector Placements:
    
    physDetector2 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta), 0, r_i, 0, 0, r_i),
        logicDetector, "physDetector", logicWorld, false, 2, true);
    
    physDetector3 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(1 * 72 * degree), pow(-1, 2) * r_i * tan(Beta) * sin(1 * 72 * degree), r_i, 0, 0, r_i ),
        logicDetector, "physDetector", logicWorld, false, 3, true);
    
    physDetector4 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), pow(-1, 2) * r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, r_i ),
        logicDetector, "physDetector", logicWorld, false, 4, true);

    physDetector5 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(1 * 72 * degree), pow(-1, 1) * r_i * tan(Beta) * sin(1* 72 * degree), r_i, 0, 0, r_i),
        logicDetector, "physDetector", logicWorld, false, 5, true);

    physDetector6 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), pow(-1, 1) * r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, r_i),
        logicDetector, "physDetector", logicWorld, false, 6, true);

    physDetector7 = new G4PVPlacement(
        doubleRotation(0, -0 * r_i * tan(Beta), 0, r_i, 0, 0, r_i ),
        logicDetector, "physDetector", logicWorld, false, 7, true);

    physDetector8 = new G4PVPlacement(
        doubleRotation(1 * 180 * degree, -1 * r_i * tan(Beta), 0, r_i, 0, 0, r_i ),
        logicDetector, "physDetector", logicWorld, false, 8, true);

    physDetector9 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(1 * 72 * degree), pow(-1, 2) * r_i * tan(Beta) * sin(1 * 72 * degree), r_i, 0, 0, r_i ),
        logicDetector, "physDetector", logicWorld, false, 9, true);

    physDetector10 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), pow(-1, 2) * r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, r_i ),
        logicDetector, "physDetector", logicWorld, false, 10, true);

    physDetector11 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(1 * 72 * degree), pow(-1, 1) * r_i * tan(Beta) * sin(1* 72 * degree), r_i, 0, 0, r_i ),
        logicDetector, "physDetector", logicWorld, false, 11, true);
    
    physDetector12 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), pow(-1, 1) * r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, r_i),
        logicDetector, "physDetector", logicWorld, false, 12, true);
    
    //kill surf
    physPanel= new G4PVPlacement(0, G4ThreeVector(0., 0., 1.5*m), logicPanel, "physPanel", logicWorld, false, 66, true);
    defineBoundaries();
    fScoringVolume = logicDoDi;
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

    G4LogicalBorderSurface *SecondSurface1 = new G4LogicalBorderSurface("SecondSurface", physSurf1, physWorld , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface2 = new G4LogicalBorderSurface("SecondSurface", physSurf2, physWorld , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface3 = new G4LogicalBorderSurface("SecondSurface", physSurf3, physWorld , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface4 = new G4LogicalBorderSurface("SecondSurface", physSurf4, physWorld , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface5 = new G4LogicalBorderSurface("SecondSurface", physSurf5, physWorld , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface6 = new G4LogicalBorderSurface("SecondSurface", physSurf6, physWorld , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface7 = new G4LogicalBorderSurface("SecondSurface", physSurf7, physWorld , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface8 = new G4LogicalBorderSurface("SecondSurface", physSurf8, physWorld , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface9 = new G4LogicalBorderSurface("SecondSurface", physSurf9, physWorld , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface10 = new G4LogicalBorderSurface("SecondSurface", physSurf10, physWorld , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface11 = new G4LogicalBorderSurface("SecondSurface", physSurf11, physWorld , SecondBoundary);
    G4LogicalBorderSurface *SecondSurface12 = new G4LogicalBorderSurface("SecondSurface", physSurf12, physWorld , SecondBoundary);
}

void MyDetectorConstruction::ConstructSDandField()
{
    MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");

    logicDetector->SetSensitiveDetector(sensDet);
}
