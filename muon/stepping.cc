#include "stepping.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction)
{
    fEventAction = eventAction;
}

MySteppingAction::~MySteppingAction()
{}

void MySteppingAction::UserSteppingAction(const G4Step *step)
{
    G4Track *track = step->GetTrack();
    G4int TrackID = track->GetTrackID();
    G4int EventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();
    
    const G4ParticleDefinition *particle = track->GetParticleDefinition();
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable(); // Load the particle table with all stored properties
    G4String particleNameMu ="mu-";
    G4String particleNamePhoton ="opticalphoton";
    G4String particleNameElectron ="e-";
    G4String particleNamePositron ="e+";
    G4ParticleDefinition *requiredMu = particleTable->FindParticle(particleNameMu);
    G4ParticleDefinition *requiredPhoton = particleTable->FindParticle(particleNamePhoton);
    G4ParticleDefinition *requiredElectron = particleTable->FindParticle(particleNameElectron);
    G4ParticleDefinition *requiredPositron = particleTable->FindParticle(particleNamePositron);
    G4AnalysisManager *man = G4AnalysisManager::Instance();

    if(particle == requiredMu)
        {
            if(!(step->GetPostStepPoint()->GetTouchable()->GetVolume() && step->GetPreStepPoint()->GetTouchable()->GetVolume()))
                return;
            G4LogicalVolume *PreVolume = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume();
            G4LogicalVolume *PostVolume = step->GetPostStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume();
            if(PreVolume->GetName() == "logicPanel" || PostVolume->GetName() == "logicDoDi")
                {
                    fEventAction->IsPanel(true);
                }
            if(PreVolume->GetName() == "logicWorld" && PostVolume->GetName() == "logicSurf")
                {            
                    G4ThreeVector stppos1 = step->GetPostStepPoint()->GetPosition();
                    fEventAction->Setpos1(stppos1);   
                    //G4cout << "Setting position 1" << G4endl;      
                }
            if(PreVolume->GetName() == "logicSurf" && PostVolume->GetName() == "logicWorld")
                {            
                    G4ThreeVector stppos2 = step->GetPostStepPoint()->GetPosition(); 
                    fEventAction->Setpos2(stppos2);  
                    //G4cout << "Setting position 2" << G4endl;                      
                }
        }
    else if(particle==requiredPhoton)
    {

        fEventAction->CountPhoton(TrackID);
           
    }
    /*else if(particle==requiredElectron)
    {
       track->SetTrackStatus(fStopAndKill);         
    }
    else if(particle==requiredPositron)
    {
       track->SetTrackStatus(fStopAndKill);         
    }*/
    //G4ThreeVector *parentID = step->GetPreStepPoint()->GetParentID();

}