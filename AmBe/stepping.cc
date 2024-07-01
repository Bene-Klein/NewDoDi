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
    const G4StepPoint* endPoint = step->GetPostStepPoint();
    const G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4int EventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();
    
    const G4ParticleDefinition *particle = track->GetParticleDefinition();
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable(); // Load the particle table with all stored properties
    G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    const MyDetectorConstruction *detectorConstruction = static_cast<const MyDetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();
    G4VProcess* process = const_cast<G4VProcess*>(endPoint->GetProcessDefinedStep());
    G4double Time = step->GetPostStepPoint()->GetGlobalTime();
    
    //
    G4String particleNameMu ="mu-";
    G4String particleNamePhoton ="opticalphoton";
    G4String particleNameElectron ="e-";
    G4String particleNamePositron ="e+";
    G4String particleNameNeutron ="neutron";
    //
    G4ParticleDefinition *requiredMu = particleTable->FindParticle(particleNameMu);
    G4ParticleDefinition *requiredPhoton = particleTable->FindParticle(particleNamePhoton);
    G4ParticleDefinition *requiredElectron = particleTable->FindParticle(particleNameElectron);
    G4ParticleDefinition *requiredPositron = particleTable->FindParticle(particleNamePositron);
    G4ParticleDefinition *requiredNeutron = particleTable->FindParticle(particleNameNeutron);

    G4AnalysisManager *man = G4AnalysisManager::Instance();    
if(!(step->GetPostStepPoint()->GetTouchable()->GetVolume() && step->GetPreStepPoint()->GetTouchable()->GetVolume()))
    return;
    G4LogicalVolume *PreVolume = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume();
    G4LogicalVolume *PostVolume = step->GetPostStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume();
if(PreVolume->GetName() == "logicPanel" || PostVolume->GetName() == "logicDoDi")
{
    fEventAction->IsPanel(true);
}
if(PreVolume->GetName() == "logicFoil" && PostVolume->GetName() == "logicWater")
{            
    G4ThreeVector stppos1 = step->GetPostStepPoint()->GetPosition();
    fEventAction->Setpos1(stppos1);       
}
if(PreVolume->GetName() == "logicWater" && PostVolume->GetName() == "logicFoil")
{            
    G4ThreeVector stppos2 = step->GetPostStepPoint()->GetPosition(); 
    fEventAction->Setpos2(stppos2);  
                        
}
if(process->GetProcessName() == "nCapture" && PreVolume->GetName() == "logicWater")
{
    fEventAction->SetProcess(1);
    fEventAction->SetTime(Time);
    //G4cout << "process:" << process->GetProcessName() << G4endl;
}
if(PreVolume->GetName() == "logicWater")
{
    G4double edep = step->GetTotalEnergyDeposit();
    fEventAction->AddEdep(edep);
}
if(particle == requiredNeutron && step->GetTrack()->GetTrackStatus() == fStopAndKill && process->GetProcessName() == "nCapture" )
    {
        G4HadronicProcess* hproc = dynamic_cast<G4HadronicProcess*>(process);
        const G4Isotope* target = NULL;
        if (hproc) target = hproc->GetTargetIsotope();
        G4String targetName = "XXXX";  
        if (target) {
            targetName = target->GetName();
            //G4cout << targetName << G4endl;
            man->FillNtupleSColumn(2,0,targetName);
            G4ThreeVector position = endPoint->GetPosition();
            G4double x_position = position.getX();
            G4double y_position = position.getY();
            G4double z_position = position.getZ();
            man->FillNtupleDColumn(2,1,x_position);
            man->FillNtupleDColumn(2,2,y_position);
            man->FillNtupleDColumn(2,3,z_position);
            man->FillNtupleIColumn(2,4,EventID);
            man->AddNtupleRow(2);  

            G4double energy = preStepPoint->GetKineticEnergy();
            //man->FillNtupleDColumn(2,4,energy);
            if(targetName.contains("Gd"))
            {
                fEventAction->Addn(1);
                fEventAction->SetTime(Time);
            }
            if((targetName.contains("H1") || targetName.contains("O16")) && PostVolume->GetName() == "logicGadolinium")
            {
                fEventAction->Addn(2);
                fEventAction->SetTime(Time);
            }
            if((targetName.contains("H1") || targetName.contains("O16")) && PostVolume->GetName() != "logicGadolinium")            
            {
                fEventAction->Addn(3);
                fEventAction->SetTime(Time);
            }
            if(targetName.contains("Mn") ||targetName.contains("Si") ||targetName.contains("Cr") ||
                targetName.contains("Ni") ||targetName.contains("Fe"))
            {
                fEventAction->Addn(4);
                fEventAction->SetTime(Time);
            }
        }
    }
}