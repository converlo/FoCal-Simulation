// This code is part of the ALICE FOCAL project, which is a forward calorimeter for the ALICE experiment at CERN. The code is written in C++ and uses the ROOT framework for data analysis and the FastJet library for jet reconstruction. The purpose of this code is to analyze jets in a grid of events, perform jet matching between reconstructed jets and truth-level jets, and store the results in a ROOT file for further analysis.




// dépendances, gestion de fichiers, données, etc
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TH3.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TROOT.h"
#include "Riostream.h"


// bibliothèques (ALICE, MC, AliFOCALCluster, FASTJET)
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliFOCALCluster.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#ifndef __FJCORE__
#include "fastjet/GhostedAreaSpec.hh" // for area support
#endif								  // __FJCORE__

using std::cout;
using std::endl;

// Forward declarations
Float_t eta(TParticle *); // Float_T = float, TParticle is a class from ROOT representing a particle in the event generator stack, and this function will calculate the pseudorapidity (eta) of the particle based on its momentum components.
Float_t phi(TParticle *); // This function will calculate the azimuthal angle (phi) of the particle based on its momentum components.

void GetMomentum(TLorentzVector &p, const Float_t *vertex, AliFOCALCluster *foClust, Float_t mass); 
// This function will calculate the momentum of a cluster (foClust) and store it in a TLorentzVector (p). It takes into account the position of the cluster, the vertex position, and an optional mass parameter.
// & = reference, * = pointer, const = constant (cannot be modified)
// TLorentzVector is a class from ROOT representing a 4-momentum, AliFOCALCluster is a class representing a cluster in the FOCAL detector.

void GetGeometricalMatchingLevel(fastjet::PseudoJet& jet1, fastjet::PseudoJet& jet2, Double_t &d);
// This function will calculate the geometrical distance (d) between two jets (jet1 and jet2) in the eta-phi space.
// fastjet::PseudoJet is a class from the FastJet library representing a jet, and Double_t = double.

void DoJetMatching(std::vector<fastjet::PseudoJet> &jetArray1, std::vector<fastjet::PseudoJet> &jetArray2);
// This function will perform jet matching between two arrays of jets (jetArray1 and jetArray2) based on their distance in the eta-phi space.
// std::vector is a C++ standard library container representing a dynamic array

class JetMatchingParams : public fastjet::PseudoJet::UserInfoBase // This class will store the parameters related to jet matching, such as the indices of the matched jets and their distances. It inherits from fastjet::PseudoJet::UserInfoBase, which allows us to attach this information to each jet in the FastJet framework. public means that the members of the class are accessible from outside the class.
{
public: // constructeur, public members can be accessed from outside the class
	JetMatchingParams() = default;
	JetMatchingParams(const int &index1, const int &index2, const double &distance1, const double &distance2) : mIndex1(index1), mIndex2(index2), mDistance1(distance1), mDistance2(distance2) {}
	JetMatchingParams(JetMatchingParams &params) : mIndex1(params.index1()), mIndex2(params.index2()), mDistance1(params.distance1()), mDistance2(params.distance2()) {}
	JetMatchingParams(const JetMatchingParams &params) : mIndex1(params.index1()), mIndex2(params.index2()), mDistance1(params.distance1()), mDistance2(params.distance2()) {}

	/// access to the jet index
	int index1() const { return mIndex1; }
	int index2() const { return mIndex2; }

	/// access to the distance
	double distance1() const { return mDistance1; }
	double distance2() const { return mDistance2; }


// pourquoi y'a deux distances une associé a chcun des deux jets je pensais qu'on rgzrdait la distance entre les deux jets pour faire le matching.


	/// Set the indices and distance
	void setJetIndexDistance1(int index, double distance)
	{
		mIndex1 = index;
		mDistance1 = distance;
	}
	void setJetIndexDistance2(int index, double distance)
	{
		mIndex2 = index;
		mDistance2 = distance;
	}

protected:
	int mIndex1 = -1;		 // jet1 index
	int mIndex2 = -1;		 // jet2 index
	double mDistance1 = 999; // Distance to jet1
	double mDistance2 = 999; // Distance to jet2
};

class SW_JetArea : public fastjet::SelectorWorker // This class defines a custom selector for jets based on their area. It inherits from fastjet::SelectorWorker, which is a base class for defining custom selectors in the FastJet framework. The selector will pass jets that have an area greater than a specified minimum area (mMinArea).
{

// area = surfcae occupé dans l'espace de phase eta, phi par le jet, plus le jet est grand plus son area est grande, on peut faire du background subtraction en soustrayant rho*area du pt du jet, ou en faisant du jet grooming en supprimant les sous-jets qui ont une area trop petite
public:
	SW_JetArea(const double &minArea) : mMinArea(minArea) {}

	// the selector's description
	std::string description() const
	{
		std::ostringstream oss;
		oss << "Minimum area " << mMinArea;
		return oss.str();
	}

	bool pass(const fastjet::PseudoJet &p) const
	{
		return p.area() > mMinArea;
	}

private:
	double mMinArea; // the vertex number we're selecting
};

fastjet::Selector SelectorMinArea(const double &minArea)
{
	return fastjet::Selector(new SW_JetArea(minArea));
}

void AnalyzeJetsGrid()
{
	const Float_t Ecut = 2.0;
	const Float_t kPi0Mass = 0.135;
	const Float_t kPiPlusMass = 0.139;

	const Int_t nR = 3;
	const Float_t Rvals[nR] = {0.2, 0.4, 0.6}; // Cone radii
	// const Float_t Rvals[nR] = {0.2, 0.3, 0.4, 0.5, 0.6}; // Cone radii

	if (gSystem->AccessPathName("focalClusters.root"))
	{
		cout << "AnalyzeJetsGrid() ERROR: focalClusters.root not found!" << endl;
		return;
	}
	TFile *clusterFile = new TFile("focalClusters.root");
	cout << "Clusters from: " << clusterFile->GetName() << endl;

	TFile *fout = new TFile("analysisJets.root", "RECREATE");
	fout->cd();

	Float_t jetE, jetpT, jetpT_sub, jetPhi, jetEta, jetRap;
	Float_t jetE_match, jetpT_match, jetpT_sub_match, jetPhi_match, jetEta_match, jetRap_match;
	Float_t jet_distmatch;
	Int_t jetR;
	Int_t ievt;

	TTree *results = new TTree("jetTree", "jets Info Tree");
	results->Branch("ievt", &ievt, "ievt/I");
	results->Branch("jetE", &jetE, "jetE/F");
	results->Branch("jetpT", &jetpT, "jetpT/F");
	results->Branch("jetpT_sub", &jetpT_sub, "jetpT_sub/F");
	results->Branch("jetPhi", &jetPhi, "jetPhi/F");
	results->Branch("jetEta", &jetEta, "jetEta/F");
	results->Branch("jetRap", &jetRap, "jetRap/F");
	results->Branch("jetE_match", &jetE_match, "jetE_match/F");
	results->Branch("jetpT_match", &jetpT_match, "jetpT_match/F");
	results->Branch("jetPhi_match", &jetPhi_match, "jetPhi_match/F");
	results->Branch("jetEta_match", &jetEta_match, "jetEta_match/F");
	results->Branch("jetRap_match", &jetRap_match, "jetRap_match/F");
	results->Branch("jet_distmatch", &jet_distmatch, "jet_distmatch/F");
	results->Branch("jetR", &jetR, "jetR/I");

	Float_t truthjetE, truthjetpT, truthjetPhi, truthjetEta, truthjetRap;
	Int_t truthjetR;

	TTree *results_truth = new TTree("truthjetTree", "MC jets Info Tree");
	results_truth->Branch("truthjetE", &truthjetE, "truthjetE/F");
	results_truth->Branch("truthjetpT", &truthjetpT, "truthjetpT/F");
	results_truth->Branch("truthjetPhi", &truthjetPhi, "truthjetPhi/F");
	results_truth->Branch("truthjetEta", &truthjetEta, "truthjetEta/F");
	results_truth->Branch("truthjetRap", &truthjetRap, "truthjetRap/F");
	results_truth->Branch("truthjetR", &truthjetR, "truthjetR/I");

	TH1D *nevt_h = new TH1D("nevt_h", "Number of Events", 1, 0, 1);

	// Alice run loader
	AliRunLoader *fRunLoader = AliRunLoader::Open("root_archive.zip#galice.root");

	if (!fRunLoader->GetAliRun())
		fRunLoader->LoadgAlice();
	if (!fRunLoader->TreeE())
		fRunLoader->LoadHeader();
	if (!fRunLoader->TreeK())
		fRunLoader->LoadKinematics();

	TObjArray primTracks;

	std::vector<fastjet::JetDefinition> jet_def;
	std::vector<fastjet::JetDefinition> jet_defBkg;
	Double_t ghost_maxrap = 6.0;
	fastjet::GhostedAreaSpec area_spec(ghost_maxrap);
	fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

	for (Int_t iR = 0; iR < nR; iR++)
	{
		jet_def.push_back(fastjet::JetDefinition(fastjet::antikt_algorithm, Rvals[iR]));
		jet_defBkg.push_back(fastjet::JetDefinition(fastjet::kt_algorithm, Rvals[iR]));
	}

	for (ievt = 0; ievt < fRunLoader->GetNumberOfEvents(); ievt++)
	{

		// Get MC Stack
		Int_t ie = fRunLoader->GetEvent(ievt);

		if (ie != 0)
			continue;
		nevt_h->Fill(0.5); // count events

		AliStack *stack = fRunLoader->Stack();

		float Vertex[3] = {0x0};

		std::vector<fastjet::PseudoJet> input_Particles;

		for (Int_t iTrk = 0; iTrk < stack->GetNtrack(); iTrk++)
		{
			TParticle *part = stack->Particle(iTrk);

			if (iTrk == 6)
			{
				Vertex[0] = part->Vx();
				Vertex[1] = part->Vy();
				Vertex[2] = part->Vz();
			}

			if (!stack->IsPhysicalPrimary(iTrk))
				continue;
			if (part->GetFirstDaughter() >= 0 && stack->IsPhysicalPrimary(part->GetFirstDaughter())) // if decay daughters are phys prim, only use those
				continue;

			if (part->Eta() > 5.8 || part->Eta() < 3.2)
				continue;

			input_Particles.push_back(fastjet::PseudoJet(part->Px(), part->Py(), part->Pz(), part->Energy()));
		}

		// Get Clusters
		TTree *tClusters = 0;
		if (clusterFile->GetDirectory(Form("Event%i", ievt)))
		{
			clusterFile->GetDirectory(Form("Event%i", ievt))->GetObject("fTreeR", tClusters);
		}
		else
		{
			cout << "Cannot find event " << ievt << " in cluster file " << clusterFile->GetName() << endl;
			clusterFile->ls();
		}

		TBranch *bClusters;
		bClusters = tClusters->GetBranch("AliFOCALCluster"); // Branch for final ECAL clusters

		TClonesArray *clustersArray = 0;
		bClusters->SetAddress(&clustersArray);
		bClusters->GetEvent(0);

		std::vector<fastjet::PseudoJet> input_Clusters;

		for (int iClust = 0; iClust < clustersArray->GetEntries(); iClust++)
		{
			AliFOCALCluster *clust = (AliFOCALCluster *)clustersArray->At(iClust);

			TLorentzVector mom;
			Float_t vertex[3] = {0.};
			GetMomentum(mom, vertex, clust, 0.);
			input_Clusters.push_back(fastjet::PseudoJet(mom.Px(), mom.Py(), mom.Pz(), clust->E()));
		}

		TBranch *AllClusters = tClusters->GetBranch("AliFOCALClusterItr");

		TClonesArray *FullclustersArray = 0;
		AllClusters->SetAddress(&FullclustersArray);
		AllClusters->GetEvent(0);

		// Selecting HCAL clusters
		for (int iClust = 0; iClust < FullclustersArray->GetEntries(); iClust++)
		{
			AliFOCALCluster *clust = (AliFOCALCluster *)FullclustersArray->At(iClust);
			if (clust->Segment() != 6)
			{
				continue;
			}

			if (clust->E() > 6500)
			{
				std::cout << "Rejecting fake cluster with large energy\n";
				continue;
			}

			TLorentzVector mom;
			Float_t vertex[3] = {0.};
			GetMomentum(mom, vertex, clust, kPiPlusMass);
			input_Clusters.push_back(fastjet::PseudoJet(mom.Px(), mom.Py(), mom.Pz(), clust->E()));
		}

		for (Int_t iR = 0; iR < nR; iR++)
		{

			// These cuts are not used for now
			Double_t AreaCut = 0.6 * TMath::Pi() * TMath::Power(Rvals[iR], 2);
			fastjet::Selector selJet = SelectorMinArea(AreaCut) && fastjet::SelectorEMin(Ecut) && fastjet::SelectorEtaRange(3.5 + Rvals[iR], 5.5 - Rvals[iR]);

			fastjet::ClusterSequenceArea clust_seq(input_Clusters, jet_def[iR], area_def);
			fastjet::ClusterSequenceArea clust_seq_truth(input_Particles, jet_def[iR], area_def);
			double ptmin = 1.0; // GeV/c, might be too large in forward regions
			std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
			std::vector<fastjet::PseudoJet> inclusive_jets_truth = sorted_by_pt(clust_seq_truth.inclusive_jets(ptmin));

			fastjet::ClusterSequenceArea clust_seqBkg(input_Clusters, jet_defBkg[iR], area_def);
			std::vector<fastjet::PseudoJet> inclusive_jetsBkg = sorted_by_pt(clust_seqBkg.inclusive_jets(0.));

			if (inclusive_jets.size() == 0)
				continue;

			Double_t Rho = 0.;

			// Perp Cone for underlying Event Subtraction
			Double_t PerpConePhi = inclusive_jets[0].phi() + TMath::Pi() / 2;
			PerpConePhi = (PerpConePhi > 2 * TMath::Pi()) ? PerpConePhi - 2 * TMath::Pi() : PerpConePhi; // fit to 0 < phi < 2pi

			Double_t PerpConeEta = inclusive_jets[0].eta();
			Double_t PerpConePt(0);

			for (unsigned int j = 0; j < input_Clusters.size(); j++)
			{

				Double_t deltaR(0);

				Double_t dPhi = TMath::Abs(input_Clusters[j].phi() - PerpConePhi);
				dPhi = (dPhi > TMath::Pi()) ? 2 * TMath::Pi() - dPhi : dPhi;

				Double_t dEta = TMath::Abs(input_Clusters[j].eta() - PerpConeEta);

				deltaR = TMath::Sqrt(TMath::Power(dEta, 2) + TMath::Power(dPhi, 2));

				if (deltaR <= Rvals[iR])
					PerpConePt += input_Clusters[j].pt();
			}

			Double_t PerpConeRho = PerpConePt / (TMath::Pi() * TMath::Power(Rvals[iR], 2));
			Rho = PerpConeRho;

			for (unsigned int i = 0; i < inclusive_jets_truth.size(); i++)
			{
				inclusive_jets_truth[i].set_user_index(i);
				inclusive_jets_truth[i].set_user_info(new JetMatchingParams());
				truthjetE = inclusive_jets_truth[i].E();
				truthjetpT = inclusive_jets_truth[i].perp();
				truthjetPhi = inclusive_jets_truth[i].phi();
				truthjetEta = inclusive_jets_truth[i].eta();
				truthjetRap = inclusive_jets_truth[i].rap();
				truthjetR = int(Rvals[iR] * 10);

				results_truth->Fill();
			}

			for (unsigned int i = 0; i < inclusive_jets.size(); i++)
			{
				inclusive_jets[i].set_user_index(i);
				inclusive_jets[i].set_user_info(new JetMatchingParams());
			}

			DoJetMatching(inclusive_jets, inclusive_jets_truth);

			for (unsigned int i = 0; i < inclusive_jets.size(); i++)
			{
				jetE = inclusive_jets[i].E();
				jetpT = inclusive_jets[i].perp();
				jetpT_sub = inclusive_jets[i].perp() - (Rho * inclusive_jets[i].area());
				jetPhi = inclusive_jets[i].phi();
				jetEta = inclusive_jets[i].eta();
				jetRap = inclusive_jets[i].rap();
				jetR = int(Rvals[iR] * 10);

				fastjet::PseudoJet *matchedJet;
				if (inclusive_jets[i].user_info<JetMatchingParams>().index1() >= 0)
				{
					matchedJet = &(inclusive_jets_truth[inclusive_jets[i].user_info<JetMatchingParams>().index1()]);
				}
				else
				{
					// Assign an empty jets in case of no match
					matchedJet = new fastjet::PseudoJet();
				}

				jetE_match = matchedJet->E();
				jetpT_match = matchedJet->perp();
				jetPhi_match = matchedJet->phi();
				jetEta_match = matchedJet->eta();
				jetRap_match = matchedJet->rap();
				jet_distmatch = inclusive_jets[i].user_info<JetMatchingParams>().distance1();

				results->Fill();
			}
		}
	}

	fout->Write();
	fout->Close();
	clusterFile->Close();
	fRunLoader->Delete();
}

Float_t eta(TParticle *part)
{
	Double_t pt = sqrt(part->Px() * part->Px() + part->Py() * part->Py());
	if (pt == 0)
		return 9999;
	return -log(tan(TMath::ATan2(pt, part->Pz()) / 2));
}

Float_t phi(TParticle *part)
{
	return TMath::ATan2(part->Py(), -part->Px());
}

//_______________________________________________________________________
void GetMomentum(TLorentzVector &mom, const Float_t *vertex, AliFOCALCluster *foClust, Float_t mass)
{
	if (mass < 0)
		mass = 0;

	Double32_t energy = foClust->E();

	if (energy < mass)
		mass = 0;

	Double_t p = TMath::Sqrt(energy * energy - mass * mass);

	// std::cout<<"Total energy: "<<energy<<" HCAL Energy is: "<<foClust->GetSegmentEnergy(0)<<" "<<foClust->GetSegmentEnergy(1)<<" "<<foClust->GetSegmentEnergy(2)<<" "<<foClust->GetSegmentEnergy(3)<<" "<<foClust->GetSegmentEnergy(4)<<" "<<foClust->GetSegmentEnergy(5)<<" "<<foClust->GetSegmentEnergy(6)<<" "<< foClust->GetHCALEnergy()<<std::endl;
	Float_t pos[3] = {foClust->X(), foClust->Y(), foClust->Z()};

	if (vertex)
	{ // calculate direction from vertex
		pos[0] -= vertex[0];
		pos[1] -= vertex[1];
		pos[2] -= vertex[2];
	}

	Double_t r = TMath::Sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);

	mom.SetPxPyPzE(p * pos[0] / r, p * pos[1] / r, p * pos[2] / r, energy);
}

void DoJetMatching(std::vector<fastjet::PseudoJet> &jetArray1, std::vector<fastjet::PseudoJet> &jetArray2)
{
	for (auto &jet1 : jetArray1)
	{
		JetMatchingParams matchingParams1(jet1.user_info<JetMatchingParams>());

		for (auto &jet2 : jetArray2)
		{
			JetMatchingParams matchingParams2(jet2.user_info<JetMatchingParams>());

			// This calculates the distance in the rap - phi 
			double dist = jet1.delta_R(jet2);
			// In order to do calculation in eta-Phi uncomment this line
			// GetGeometricalMatchingLevel(jet1, jet2, dist);

			if (dist < matchingParams1.distance1())
			{
				matchingParams1.setJetIndexDistance2(matchingParams1.index1(), matchingParams1.distance1());
				matchingParams1.setJetIndexDistance1(jet2.user_index(), dist);
			}
			else if (dist < matchingParams1.distance2())
			{
				matchingParams1.setJetIndexDistance2(jet2.user_index(), dist);
			}

			if (dist < matchingParams2.distance1())
			{
				matchingParams2.setJetIndexDistance2(matchingParams2.index1(), matchingParams2.distance1());
				matchingParams2.setJetIndexDistance1(jet1.user_index(), dist);
			}
			else if (dist < matchingParams2.distance2())
			{
				matchingParams2.setJetIndexDistance2(jet1.user_index(), dist);
			}

			jet2.user_info_shared_ptr() = fastjet::SharedPtr<fastjet::PseudoJet::UserInfoBase>(new JetMatchingParams(matchingParams2));
		}

		jet1.user_info_shared_ptr() = fastjet::SharedPtr<fastjet::PseudoJet::UserInfoBase>(new JetMatchingParams(matchingParams1));
	}

	for (auto &jet1 : jetArray1)
	{
		fastjet::PseudoJet *jet2 = jet1.user_info<JetMatchingParams>().index1() >= 0 ? &(jetArray2[jet1.user_info<JetMatchingParams>().index1()]) : nullptr;
		fastjet::PseudoJet *jet3 = jet1.user_info<JetMatchingParams>().index2() >= 0 ? &(jetArray2[jet1.user_info<JetMatchingParams>().index2()]) : nullptr;

		if (jet2 && jet2->user_info<JetMatchingParams>().index1() != jet1.user_index())
		{
			if (!jet3 || jet3->user_info<JetMatchingParams>().index1() != jet1.user_index())
			{
				JetMatchingParams params;
				jet1.user_info_shared_ptr() = fastjet::SharedPtr<fastjet::PseudoJet::UserInfoBase>(new JetMatchingParams(params));
			}
			else
			{
				JetMatchingParams matchingParams(jet1.user_info<JetMatchingParams>());
				matchingParams.setJetIndexDistance1(jet1.user_info<JetMatchingParams>().index2(), jet1.user_info<JetMatchingParams>().distance2());
				matchingParams.setJetIndexDistance2(-1, 999);
				jet1.user_info_shared_ptr() = fastjet::SharedPtr<fastjet::PseudoJet::UserInfoBase>(new JetMatchingParams(matchingParams));
			}
		}
	}
}

void GetGeometricalMatchingLevel(fastjet::PseudoJet& jet1, fastjet::PseudoJet& jet2, Double_t &d)
{
	Double_t deta = jet2.eta() - jet1.eta();
	Double_t dphi = jet2.phi() - jet1.phi();
	dphi = TVector2::Phi_mpi_pi(dphi);
	d = TMath::Sqrt(deta * deta + dphi * dphi);
}