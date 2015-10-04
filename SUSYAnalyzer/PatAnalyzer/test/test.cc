MiniAODGenPartAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
        using namespace edm;
        using namespace reco;
        using namespace pat;

        // Pruned particles are the one containing "important" stuff
        Handle<edm::View<reco::GenParticle> > pruned;
        iEvent.getByToken(prunedGenToken_,pruned);

        // Packed particles are all the status 1, so usable to remake jets
        // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
        Handle<edm::View<pat::PackedGenParticle> > packed;
        iEvent.getByToken(packedGenToken_,packed);

        //let's try to find all status1 originating directly from a B meson decay 

        for(size_t i=0; i<pruned->size();i++){
                if(abs((*pruned)[i].pdgId()) > 500 && abs((*pruned)[i].pdgId()) <600){
                        const Candidate * bMeson = &(*pruned)[i];
                        std::cout << "PdgID: " << bMeson->pdgId() << " pt " << bMeson->pt() << " eta: " << bMeson->eta() << " phi: " << bMeson->phi() << std::endl;
                        std::cout << "  found daugthers: " << std::endl;
                        for(size_t j=0; j<packed->size();j++){
//get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection 
                                const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
                                if(motherInPrunedCollection != nullptr && isAncestor( bMeson , motherInPrunedCollection)){
                                        std::cout << "     PdgID: " << (*packed)[j].pdgId() << " pt " << (*packed)[j].pt() << " eta: " << (*packed)[j].eta() << " phi: " << (*packed)[j].phi() << std::endl;
                                }
                        }
                }

        }


}

bool MiniAODGenPartAnalyzer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
//particle is already the ancestor
        if(ancestor == particle ) return true;

//otherwise loop on mothers, if any and return true if the ancestor is found
        for(size_t i=0;i< particle->numberOfMothers();i++)
        {
                if(isAncestor(ancestor,particle->mother(i))) return true;
        }
//if we did not return yet, then particle and ancestor are not relatives
        return false;
}
