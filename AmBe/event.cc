#include "event.hh"


MyEventAction::MyEventAction(MyRunAction*)
{
    fEdep = false;   

    length = 0.;
}

MyEventAction::~MyEventAction()
{}

void MyEventAction::BeginOfEventAction(const G4Event*)
{
    fEdep = false;   
    fpos1 = G4ThreeVector(0.,0.,0.);
    fpos2 = G4ThreeVector(0.,0.,0.);
    isPanel = false;
    length = 0.;
    counter = 0;
    Tracks.clear();
    prnumber = 0;
    process = 0;    
    CaptureTime =0;
    n =0;
    nWat =0;
    nSteel =0;
}

    


void MyEventAction::EndOfEventAction(const G4Event* event)
{

    G4int EventID = event->GetEventID();

    if (EventID % 100 == 0) {
            end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end - start;
            start=end;
            std::time_t end_time = std::chrono::system_clock::to_time_t(end);
            #ifndef G4MULTITHREADED
            G4cout << "ID: "<< EventID << "  Time: " << elapsed_seconds.count()  << "  sysTime: " << std::ctime(&end_time)<< G4endl;
            #endif
        }
    if (EventID % 100000 == 0) 
    {
        G4cout << "ID: "<< EventID <<G4endl;
    }
    G4double X1=fpos1.x();
    G4double X2=fpos2.x();
    G4double Y1=fpos1.y();
    G4double Y2=fpos2.y();
    G4double Z1=fpos1.z();
    G4double Z2=fpos2.z(); 

        if((std::find(Tracks.begin(), Tracks.end(), totPhoton) == Tracks.end())) 
        {
            #ifndef G4MULTITHREADED
            //G4cout<< "TrackID:"<< Tracks.size() << G4endl;
            #endif
            Tracks.insert(Tracks.begin(),totPhoton);

        }


    G4double len = sqrt(pow((X1-X2),2)+pow((Y1-Y2),2)+pow((Z1-Z2),2));
    #ifndef G4MULTITHREADED
    //G4cout << "Length: "<< len << G4endl;
    //G4cout << "Energy deposition: " << fEdep << G4endl;
    #endif

    
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    

    man->FillNtupleDColumn(1,0,CaptureTime);
    man->FillNtupleIColumn(1,1,n);
    man->FillNtupleIColumn(1,2,EventID);
    man->AddNtupleRow(1);  

}
