#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4AnalysisManager.hh"
#include <chrono>
#include <ctime>
#include "run.hh"

class MyEventAction : public G4UserEventAction
{
    public:
        MyEventAction(MyRunAction*);
        ~MyEventAction();

        virtual void BeginOfEventAction(const G4Event*);
        virtual void EndOfEventAction(const G4Event*);

        void Setpos1(G4ThreeVector pos1) {fpos1 = pos1;};
        void Setpos2(G4ThreeVector pos2) {fpos2 = pos2;};
        void CountPhoton(G4int tracks) {totPhoton = tracks;};
        void IsPanel(G4bool panel) {isPanel=panel;};
        void AddEdep(G4double edep) { fEdep += edep; } 
        void SetProcessName(G4String processName) {ProcessName=processName;};
        void SetProcess(G4int processName) {process=processName;};
        void SetTime(G4double capTime) {CaptureTime=capTime;};
        void Addn(G4int i) {n=i;};    
    private:
        G4bool fEdep;
        G4int EventID;
        G4int totPhoton;
        G4int counter;
        G4int prnumber;
        G4int process;
        G4int n;
        G4int nWat;
        G4int nSteel;
        G4double length;
        G4double CaptureTime;
        G4String ProcessName;
        G4ThreeVector fpos1;
        G4ThreeVector fpos2;
        G4bool isPanel;

        std::vector<G4int> Tracks;
        std::chrono::time_point<std::chrono::system_clock> end;
        std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
};


#endif