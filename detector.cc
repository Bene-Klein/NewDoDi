#include "detector.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{
    quEff = new G4PhysicsOrderedFreeVector();

    std::ifstream datafile;
    datafile.open("eff.dat");

    while(1)
    {
        G4double wlen, queff;

        datafile >> wlen >> queff;

        if(datafile.eof())
            break;

        G4cout << wlen << " " << queff << G4endl;

        quEff-> InsertValues(wlen, queff/100.);
    }
    datafile.close();
}
MySensitiveDetector::~MySensitiveDetector()
{}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
    
    G4Track *track = aStep->GetTrack();
    const G4ParticleDefinition *particle = track->GetParticleDefinition();
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable(); // Load the particle table with all stored properties
    G4String particleNameopt ="opticalphoton";
    G4String particleNameNeut ="neutron";

    G4ParticleDefinition *requiredParticle = particleTable->FindParticle(particleNameopt);
    G4ParticleDefinition *requiredParticleNeut = particleTable->FindParticle(particleNameNeut);



    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

    G4ThreeVector posPhoton = preStepPoint->GetPosition();

    //G4cout << "positition: "<< posPhoton << G4endl;

    const G4VTouchable *touchable =aStep->GetPreStepPoint()->GetTouchable();

    G4int copyNo = touchable->GetCopyNumber();

    G4int parentId = aStep->GetTrack()->GetParentID();
    
    G4int EventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();

    G4double Time = aStep->GetPreStepPoint()->GetGlobalTime();
        
    G4VPhysicalVolume *physVol = touchable->GetVolume();
    G4ThreeVector posDetector = physVol->GetTranslation();
    G4ThreeVector momPhoton = preStepPoint->GetMomentum();

    G4double wlen = (1.239841939*eV/momPhoton.mag())*1E+03;
    //G4cout << "Detector position: "<< posDetector <<G4endl;

    G4AnalysisManager *man = G4AnalysisManager::Instance();
    //G4cout<<" Detector"<<G4endl;
    
if(particle == requiredParticle)
{
    track->SetTrackStatus(fStopAndKill); 
    if(G4UniformRand() < quEff ->Value(wlen) )
    {   
        man->FillNtupleDColumn(0,0,Time);
        man->FillNtupleIColumn(0,1,copyNo);
        man->FillNtupleIColumn(0,2,EventID); //FillNtupleDColumn(Ntuple Number, Entry Number, Entry)
        man->AddNtupleRow(0);    
    }
}
    return true;
}