#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4Tet.hh"
#include "G4Polyhedra.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
#include "G4RotationMatrix.hh"
#include "G4MultiUnion.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh" 
#include "G4GenericMessenger.hh"
#include "G4OpticalSurface.hh"
#include "G4Material.hh"
#include "G4MultiUnion.hh"
#include "G4UnionSolid.hh"
#include "G4UserLimits.hh"
#include <cmath>

#include "detector.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public: 
    MyDetectorConstruction();
    ~MyDetectorConstruction();

    G4LogicalVolume *GetScoringVolume() const { return fScoringVolume; }

    virtual G4VPhysicalVolume *Construct();
private: 
    virtual void ConstructSDandField();

    
    G4Transform3D Rotation(G4double, G4double, G4double, G4double, G4double, G4double, G4double);
    G4Transform3D doubleRotation(G4double, G4double, G4double, G4double, G4double, G4double, G4double);
    G4Transform3D rotZ(G4double theta, G4double x_1, G4double y_1, G4double z_1);
    G4LogicalVolume* MyDoDiConstruction(G4double a,G4double Steel,G4double Foil,G4Material* mat, char* name);
    void SetMaxStep(G4double maxStep);
    
    G4int nCols, nRows;
    G4double gerscale,thetax,thetay,thetaz;

    G4Polyhedra *solidWater, *solidFoil, *solidSteel,*solidWater1, *solidFoil1, *solidSteel1, *solidGadolinium;
    G4Polycone *solidDetector,*solidTube;
    G4Box *solidGermanium, *solidWorld, *solidPanel, *solidParaffin;
    G4Tubs *solidHuman;
    G4LogicalVolume *logicWorld,*logicWater,*logicWater1, *logicGermanium, *logicDetector,
                    *logicObject, *logicmesh, *logicTube, *logicPanel, *logicFoil,*logicFoil1, *logicParaffin,
                    *logicSteel,*logicSteel1, *logicHuman, *logicGadolinium;
    G4VPhysicalVolume   *physWorld, *physTube, *physPanel,*physHuman,
                        *physSteel1,*physSteel2,*physSteel3,*physSteel4,*physSteel5, *physSteel6,*physSteel7,*physSteel8,*physSteel9,*physSteel10,*physSteel11,*physSteel12,
                        *physWater1,*physWater2,*physWater3,*physWater4,*physWater5, *physWater6,*physWater7,*physWater8,*physWater9,*physWater10,*physWater11,*physWater12,
                        *physFoil1, *physFoil2, *physFoil3, *physFoil4, *physFoil5, *physFoil6, *physFoil7, *physFoil8, *physFoil9, *physFoil10, *physFoil11, *physFoil12,
                        *physDetector2, *physDetector3, *physDetector4, *physDetector5, *physDetector6, *physDetector7, *physDetector8, *physDetector9, *physDetector10, *physDetector11, *physDetector12,
                        *physParaffin1,*physParaffin2,*physParaffin3,*physParaffin4,*physParaffin5,
                        *physGadolinium;
    G4LogicalBorderSurface *waterSurface1,*waterSurface2,*waterSurface3,*waterSurface4,
    *waterSurface5,*waterSurface6,*waterSurface7,*waterSurface8,
    *waterSurface9,*waterSurface10,*waterSurface11,*waterSurface12;

    G4Material *worldMat, *germanium, *water;
    G4GenericMessenger *fMessenger;

    void DefineMaterials();
    void defineBoundaries();
    G4Material *fWLSfoilPMMA, *matsteel, *fH2O, *fworldMat, *matparaffin, *Gdsol,*Gdsol2,  *Gdmat;
    std::vector<G4double> fenergySmall,fenergy;
    G4Element *elGd,*H, *C, *O;

    G4LogicalVolume *fScoringVolume;

    //G4int nRows, nCols;
};
#endif
