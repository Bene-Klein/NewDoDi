#include "run.hh"

MyRunAction::MyRunAction()
{   
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->CreateNtuple("Time","Time");
    man->CreateNtupleDColumn("Time");
    man->CreateNtupleIColumn("Detector");
    man->CreateNtupleIColumn("EventID");
    man->FinishNtuple(0); 

    man->CreateNtuple("Length","Length");
    man->CreateNtupleDColumn("Length");
    man->CreateNtupleIColumn("totalPhoton");
    man->CreateNtupleIColumn("isPanel");
    man->CreateNtupleIColumn("EventID");
    man->FinishNtuple(1);
    /*  
    man->CreateNtuple("trackid","trackid");
    man->CreateNtupleDColumn("TrackID");
    man->CreateNtupleIColumn("EventID");
    man->FinishNtuple(2);

      
    man->CreateNtuple("Count","Count");
    man->CreateNtupleIColumn("Count");
    man->CreateNtupleIColumn("EventID");
    man->FinishNtuple(2);
    */    
}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();

    G4int runID = run->GetRunID();

    std::stringstream strRunID;
    strRunID << runID;

    man->OpenFile("output"+strRunID.str()+".csv");
}

void MyRunAction::EndOfRunAction(const G4Run*)
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();

    man->Write();
    man->CloseFile();
}