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

        
    private:
        G4bool fEdep;
        G4int EventID;
        G4int totPhoton;
        G4int counter;
        G4double length;
        G4ThreeVector fpos1;
        G4ThreeVector fpos2;
        G4bool isPanel;

        std::vector<G4int> Tracks;
        std::chrono::time_point<std::chrono::system_clock> end;
        std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
};


#endif