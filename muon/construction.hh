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
#include "CADMesh.hh"
#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh" 
#include "G4GenericMessenger.hh"
#include "G4OpticalSurface.hh"
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


    G4int nCols, nRows;
    G4double gerscale,thetax,thetay,thetaz;

    G4Polyhedra *solidDoDi, *solidSurf;
    G4Polycone *solidDetector,*solidTube;
    G4Box *solidGermanium, *solidWorld, *solidPanel;
    G4LogicalVolume *logicWorld,*logicDoDi,*logicDoDi1, *logicGermanium, *logicDetector, *logicObject, *logicmesh, *logicTube, *logicPanel, *logicSurf;
    G4VPhysicalVolume *physWorld,*physDoDi1,*physDoDi2,*physDoDi3,*physDoDi4,*physDoDi5,
    *physDoDi6,*physDoDi7,*physDoDi8,*physDoDi9,*physDoDi10,*physDoDi11,*physDoDi12,
    *physDetector2,*physDetector3,*physDetector4,*physDetector5,*physDetector6,
    *physDetector7,*physDetector8,*physDetector9,*physDetector10,*physDetector11,*physDetector12,
    *physTube, *physPanel, *physSurf1, *physSurf2, *physSurf3, *physSurf4, *physSurf5, *physSurf6, *physSurf7, *physSurf8,
    *physSurf9, *physSurf10, *physSurf11, *physSurf12;
    G4LogicalBorderSurface *waterSurface1,*waterSurface2,*waterSurface3,*waterSurface4,
    *waterSurface5,*waterSurface6,*waterSurface7,*waterSurface8,
    *waterSurface9,*waterSurface10,*waterSurface11,*waterSurface12;

    G4Material *worldMat, *germanium, *water;
    G4GenericMessenger *fMessenger;

    void DefineMaterials();
    void defineBoundaries();
    G4Material *fWLSfoilPMMA, *fH2O, *fworldMat;
    std::vector<G4double> fenergySmall,fenergy;
    G4Element *H, *C, *O;

    G4LogicalVolume *fScoringVolume;

    //G4int nRows, nCols;
};
#endif
