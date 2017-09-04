// fstext/determinize-star-inl.h

// Copyright 2009-2011  Microsoft Corporation;  Jan Silovsky
//           2015 Hainan Xu

// See ../../COPYING for clarification regarding multiple authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
// THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
// WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
// MERCHANTABLITY OR NON-INFRINGEMENT.
// See the Apache 2 License for the specific language governing permissions and
// limitations under the License.

#ifndef KALDI_FSTEXT_DETERMINIZE_STAR_INL_H_
#define KALDI_FSTEXT_DETERMINIZE_STAR_INL_H_
// Do not include this file directly.  It is included by determinize-star.h

#include "base/kaldi-error.h"

#ifdef _MSC_VER
#include <unordered_map>
using std::unordered_map;
#elif __cplusplus > 199711L || defined(__GXX_EXPERIMENTAL_CXX0X__)
#include <unordered_map>
using std::unordered_map;
#else
#include <tr1/unordered_map>
using std::tr1::unordered_map;
#endif

#include <vector>
#include <climits>

namespace fst {

// This class maps back and forth from/to integer id's to sequences of strings.
// used in determinization algorithm.

template<class Label, class StringId> class StringRepository {
  // Label and StringId are both integer types, possibly the same.
  // This is a utility that maps back and forth between a vector<Label> and StringId
  // representation of sequences of Labels.  It is to save memory, and to save compute.
  // We treat sequences of length zero and one separately, for efficiency.

 public:
  class VectorKey { // Hash function object.
   public:
    size_t operator()(const vector<Label> *vec) const {
      assert(vec != NULL);
      size_t hash = 0, factor = 1;
      for (typename vector<Label>::const_iterator it = vec->begin(); it != vec->end(); it++) {
        hash += factor*(*it); factor*=103333;  // just an arbitrary prime number.
      }
      return hash;
    }
  };
  class VectorEqual {  // Equality-operator function object.
   public:
    size_t operator()(const vector<Label> *vec1, const vector<Label> *vec2) const {
      return (*vec1 == *vec2);
    }
  };

  typedef unordered_map<const vector<Label>*, StringId, VectorKey, VectorEqual> MapType;

  StringId IdOfEmpty() { return no_symbol; }

  StringId IdOfLabel(Label l) {
    if (l>= 0 && l <= (Label) single_symbol_range) {
      return l + single_symbol_start;
    } else {
      // l is out of the allowed range so we have to treat it as a sequence of length one.  Should be v. rare.
      vector<Label> v; v.push_back(l);
      return IdOfSeqInternal(v);
    }
  }

  StringId IdOfSeq(const vector<Label> &v) {  // also works for sizes 0 and 1.
    size_t sz = v.size();
    if (sz == 0) return no_symbol;
    else if (v.size() == 1) return IdOfLabel(v[0]);
    else return IdOfSeqInternal(v);
  }

  inline bool IsEmptyString(StringId id) {
    return id == no_symbol;
  }
  void SeqOfId(StringId id, vector<Label> *v) {
    if (id == no_symbol) v->clear();
    else if (id>=single_symbol_start) {
      v->resize(1); (*v)[0] = id - single_symbol_start;
    } else {
      assert(static_cast<size_t>(id) < vec_.size());
      *v = *(vec_[id]);
    }
  }
  StringId RemovePrefix(StringId id, size_t prefix_len) {
    if (prefix_len == 0) return id;
    else {
      vector<Label> v;
      SeqOfId(id, &v);
      size_t sz = v.size();
      assert(sz >= prefix_len);
      vector<Label> v_noprefix(sz - prefix_len);
      for (size_t i = 0;i < sz-prefix_len;i++) v_noprefix[i] = v[i+prefix_len];
      return IdOfSeq(v_noprefix);
    }
  }

  StringRepository() {
    // The following are really just constants but don't want to complicate compilation so make them
    // class variables.  Due to the brokenness of <limits>, they can't be accessed as constants.
    string_end = (numeric_limits<StringId>::max() / 2) - 1;  // all hash values must be <= this.
    no_symbol = (numeric_limits<StringId>::max() / 2);  // reserved for empty sequence.
    single_symbol_start =  (numeric_limits<StringId>::max() / 2) + 1;
    single_symbol_range =  numeric_limits<StringId>::max() - single_symbol_start;
  }
  void Destroy() {
    for (typename vector<vector<Label>* >::iterator iter = vec_.begin(); iter != vec_.end(); ++iter)
      delete *iter;
    vector<vector<Label>* > tmp_vec;
    tmp_vec.swap(vec_);
    MapType tmp_map;
    tmp_map.swap(map_);
  }
  ~StringRepository() {
    Destroy();
  }

 private:
  DISALLOW_COPY_AND_ASSIGN(StringRepository);

  StringId IdOfSeqInternal(const vector<Label> &v) {
    typename MapType::iterator iter = map_.find(&v);
    if (iter != map_.end()) {
      return iter->second;
    } else {  // must add it to map.
      StringId this_id = (StringId) vec_.size();
      vector<Label> *v_new = new vector<Label> (v);
      vec_.push_back(v_new);
      map_[v_new] = this_id;
      assert(this_id < string_end);  // or we used up the labels.
      return this_id;
    }
  }

  vector<vector<Label>* > vec_;
  MapType map_;

  static const StringId string_start = (StringId) 0;  // This must not change.  It's assumed.
  StringId string_end;  // = (numeric_limits<StringId>::max() / 2) - 1; // all hash values must be <= this.
  StringId no_symbol;  // = (numeric_limits<StringId>::max() / 2); // reserved for empty sequence.
  StringId single_symbol_start;  // =  (numeric_limits<StringId>::max() / 2) + 1;
  StringId single_symbol_range;  // =  numeric_limits<StringId>::max() - single_symbol_start;
};


template<class F> class DeterminizerStar {
  typedef typename F::Arc Arc;
 public:
  // Output to Gallic acceptor (so the strings go on weights, and there is a 1-1 correspondence
  // between our states and the states in ofst.  If destroy == true, release memory as we go
  // (but we cannot output again).
  void Output(MutableFst<GallicArc<Arc> >  *ofst, bool destroy = true);

  // Output to standard FST.  We will create extra states to handle sequences of symbols
  // on the output.  If destroy == true, release memory as we go
  // (but we cannot output again).

  void  Output(MutableFst<Arc> *ofst, bool destroy = true);


  // Initializer.  After initializing the object you will typically call
  // Determinize() and then one of the Output functions.
  DeterminizerStar(const Fst<Arc> &ifst, float delta = kDelta,
                   int max_states = -1, bool allow_partial = false):
      ifst_(ifst.Copy()), delta_(delta), max_states_(max_states),
      determinized_(false), allow_partial_(allow_partial),
      is_partial_(false), equal_(delta),
      hash_(ifst.Properties(kExpanded, false) ?
              down_cast<const ExpandedFst<Arc>*,
              const Fst<Arc> >(&ifst)->NumStates()/2 + 3 : 20,
            hasher_, equal_),
      epsilon_closure_(ifst_, max_states, &repository_, delta) { }

  void Determinize(bool *debug_ptr) {
    assert(!determinized_);
    // This determinizes the input fst but leaves it in the "special format"
    // in "output_arcs_".
    InputStateId start_id = ifst_->Start();
    if (start_id == kNoStateId) { determinized_ = true; return; } // Nothing to do.
    else {  // Insert start state into hash and queue.
      Element elem;
      elem.state = start_id;
      elem.weight = Weight::One();
      elem.string = repository_.IdOfEmpty();  // Id of empty sequence.
      vector<Element> vec;
      vec.push_back(elem);
      OutputStateId cur_id = SubsetToStateId(vec);
      assert(cur_id == 0 && "Do not call Determinize twice.");
    }
    while (!Q_.empty()) {
      pair<vector<Element>*, OutputStateId> cur_pair = Q_.front();
      Q_.pop_front();
      ProcessSubset(cur_pair);
      if (debug_ptr && *debug_ptr) Debug();  // will exit.
      if (max_states_ > 0 && output_arcs_.size() > max_states_) {
        if (allow_partial_ == false) {
          std::cerr << "Determinization aborted since passed " << max_states_
                    << " states.\n";
          throw std::runtime_error("max-states reached in determinization");
        } else {
          KALDI_WARN << "Determinization terminated since passed " << max_states_
                     << " states, partial results will be generated.";
          is_partial_ = true;
          break;
        }
      }
    }
    determinized_ = true;
  }

  bool IsPartial() {
    return is_partial_;
  }

  // frees all except output_arcs_, which contains the important info
  // we need to output.
  void FreeMostMemory() {
    if (ifst_) {
      delete ifst_;
      ifst_ = NULL;
    }
    for (typename SubsetHash::iterator iter = hash_.begin();
        iter != hash_.end(); ++iter)
      delete iter->first;
    SubsetHash tmp;
    tmp.swap(hash_);
  }

  ~DeterminizerStar() {
    FreeMostMemory();
  }
 private:
  typedef typename Arc::Label Label;
  typedef typename Arc::Weight Weight;
  typedef typename Arc::StateId InputStateId;
  typedef typename Arc::StateId OutputStateId;  // same as above but distinguish states in output Fst.
  typedef typename Arc::Label StringId;  // Id type used in the StringRepository
  typedef StringRepository<Label, StringId> StringRepositoryType;


  // Element of a subset [of original states]

  struct Element {
    InputStateId state;
    StringId string;
    Weight weight;
    bool operator != (const Element &other) const  {
      return (state != other.state || string != other.string ||
              weight != other.weight);
    }
  };

  // Arcs in the format we temporarily create in this class (a representation, essentially of
  // a Gallic Fst).
  struct TempArc {
    Label ilabel;
    StringId ostring;  // Look it up in the StringRepository, it's a sequence of Labels.
    OutputStateId nextstate;  // or kNoState for final weights.
    Weight weight;
  };


  // Hashing function used in hash of subsets.
  // A subset is a pointer to vector<Element>.
  // The Elements are in sorted order on state id, and without repeated states.
  // Because the order of Elements is fixed, we can use a hashing function that is
  // order-dependent.  However the weights are not included in the hashing function--
  // we hash subsets that differ only in weight to the same key.  This is not optimal
  // in terms of the O(N) performance but typically if we have a lot of determinized
  // states that differ only in weight then the input probably was pathological in some way,
  // or even non-determinizable.
  //   We don't quantize the weights, in order to avoid inexactness in simple cases.
  // Instead we apply the delta when comparing subsets for equality, and allow a small
  // difference.

  class SubsetKey {
   public:
    size_t operator ()(const vector<Element> * subset) const {  // hashes only the state and string.
      size_t hash = 0, factor = 1;
      for (typename vector<Element>::const_iterator iter = subset->begin();
           iter != subset->end(); ++iter) {
        hash *= factor;
        hash += iter->state + 103333*iter->string;
        factor *= 23531;  // these numbers are primes.
      }
      return hash;
    }
  };

  // This is the equality operator on subsets.  It checks for exact match on state-id
  // and string, and approximate match on weights.
  class SubsetEqual {
   public:
    bool operator ()(const vector<Element> *s1,
                     const vector<Element> *s2) const {
      size_t sz = s1->size();
      assert(sz >= 0);
      if (sz != s2->size()) return false;
      typename vector<Element>::const_iterator iter1 = s1->begin(),
          iter1_end = s1->end(), iter2 = s2->begin();
      for (; iter1 < iter1_end; ++iter1, ++iter2) {
        if (iter1->state != iter2->state ||
           iter1->string != iter2->string ||
           ! ApproxEqual(iter1->weight, iter2->weight, delta_))
          return false;
      }
      return true;
    }
    float delta_;
    SubsetEqual(float delta): delta_(delta) {}
    SubsetEqual(): delta_(kDelta) {}
  };

  // Operator that says whether two Elements have the same states.
  // Used only for debug.
  class SubsetEqualStates {
   public:
    bool operator ()(const vector<Element> *s1, const vector<Element> *s2) const {
      size_t sz = s1->size();
      assert(sz>=0);
      if (sz != s2->size()) return false;
      typename vector<Element>::const_iterator iter1 = s1->begin(),
          iter1_end = s1->end(), iter2=s2->begin();
      for (; iter1 < iter1_end; ++iter1, ++iter2) {
        if (iter1->state != iter2->state) return false;
      }
      return true;
    }
  };

  // Define the hash type we use to store subsets.
  typedef unordered_map<const vector<Element>*, OutputStateId, SubsetKey, SubsetEqual> SubsetHash;

  class EpsilonClosure {
   public:
    EpsilonClosure(const Fst<Arc> *ifst, int max_states,
        StringRepository<Label, StringId> *repository, float delta):
      ifst_(ifst), max_states_(max_states), repository_(repository),
      delta_(delta) {

    }

    // This function computes epsilon closure of subset of states by following epsilon links.
    // Called by ProcessSubset.
    // Has no side effects except on the repository.
    void GetEpsilonClosure(const vector<Element> &input_subset,
                        vector<Element> *output_subset);

   private:
    struct EpsilonClosureInfo {
      EpsilonClosureInfo() {}
      EpsilonClosureInfo(const Element &e, const Weight &w, bool i) :
        element(e), weight_to_process(w), in_queue(i) {}
      // the weight in the Element struct is the total current weight
      // that has been processed already
      Element element;
      // this stores the weight that we haven't processed (propagated)
      Weight weight_to_process;
      // whether "this" struct is in the queue
      // we store the info here so that we don't have to look it up every time
      bool in_queue;
      bool operator<(const EpsilonClosureInfo &other) const {
        return this->element.state < other.element.state;
      }
    };

    // to further speed up EpsilonClosure() computation, we have 2 queues
    // the 2nd queue is used when we first iterate over the input set - 
    // if queue_2_.empty() then we directly set output_set equal to input_set
    // and return immediately
    // Since Epsilon arcs are relatively rare, this way we could efficiently
    // detect the epsilon-free case, without having to waste our computation e.g.
    // allocating the EpsilonClosureInfo structure; this also lets us do a 
    // level-by-level traversal, which could avoid some (unfortunately not all)
    // duplicate computation if epsilons form a DAG that is not a tree
    //
    // We put the queues here for better efficiency for memory allocation
    deque<typename Arc::StateId> queue_;
    vector<Element> queue_2_;

    // the following 2 structures together form our *virtual "map"*
    // basically we need a map from state_id to EpsilonClosureInfo that operates
    // in O(1) time, while still takes relatively small mem, and this does it well
    // for efficiency we don't clear id_to_index_ of its outdated information
    // As a result each time we do a look-up, we need to check
    // if (ecinfo_[id_to_index_[id]].element.state == id)
    // Yet this is still faster than using a std::map<StateId, EpsilonClosureInfo>
    vector<int> id_to_index_;
    // unlike id_to_index_, we clear the content of ecinfo_ each time we call
    // EpsilonClosure(). This needed because we need an efficient way to 
    // traverse the virtual map - it is just too costly to traverse the 
    // id_to_index_ vector. 
    vector<EpsilonClosureInfo> ecinfo_;

    // Add one element (elem) into cur_subset
    // it also adds the necessary stuff to queue_, set the correct weight
    void AddOneElement(const Element &elem, const Weight &unprocessed_weight);

    // Sub-routine that we call in EpsilonClosure()
    // It takes the current "unprocessed_weight" and propagate it to the 
    // states accessible from elem.state by an epsilon arc
    // and add the results to cur_subset.
    // save_to_queue_2 is set true when we iterate over the initial subset
    // - then we save it to queue_2 s.t. if it's empty, we directly return
    // the input set
    void ExpandOneElement(const Element &elem,
                          bool sorted,
                          const Weight &unprocessed_weight,
                          bool save_to_queue_2 = false);

    // no pointers below would take the ownership
    const Fst<Arc> *ifst_;
    int max_states_;
    StringRepository<Label, StringId> *repository_;
    float delta_;
  };


  // This function works out the final-weight of the determinized state.
  // called by ProcessSubset.
  // Has no side effects except on the variable repository_, and output_arcs_.

  void ProcessFinal(const vector<Element> &closed_subset, OutputStateId state) {
    // processes final-weights for this subset.
    bool is_final = false;
    StringId final_string = 0;  // = 0 to keep compiler happy.
    Weight final_weight = Weight::One();  // This value will never be accessed, and
    // we just set it to avoid spurious compiler warnings.  We avoid setting it
    // to Zero() because floating-point infinities can sometimes generate
    // interrupts and slow things down.
    typename vector<Element>::const_iterator iter = closed_subset.begin(),
        end = closed_subset.end();
    for (; iter != end; ++iter) {
      const Element &elem = *iter;
      Weight this_final_weight = ifst_->Final(elem.state);
      if (this_final_weight != Weight::Zero()) {
        if (!is_final) {  // first final-weight
          final_string = elem.string;
          final_weight = Times(elem.weight, this_final_weight);
          is_final = true;
        } else {  // already have one.
          if (final_string != elem.string) {
            std::cerr << "DeterminizerStar: FST was not functional -> not determinizable\n";
            throw std::runtime_error("Non-functional FST: cannot determinize.\n");
          }
          final_weight = Plus(final_weight, Times(elem.weight, this_final_weight));
        }
      }
    }
    if (is_final) {
      // store final weights in TempArc structure, just like a transition.
      TempArc temp_arc;
      temp_arc.ilabel = 0;
      temp_arc.nextstate = kNoStateId;  // special marker meaning "final weight".
      temp_arc.ostring = final_string;
      temp_arc.weight = final_weight;
      output_arcs_[state].push_back(temp_arc);
    }
  }

  // ProcessTransition is called from "ProcessTransitions".  Broken out for
  // clarity.  Has side effects on output_arcs_, and (via SubsetToStateId), Q_
  // and hash_.
  void ProcessTransition(OutputStateId state, Label ilabel, vector<Element> *subset);

  // "less than" operator for pair<Label, Element>.   Used in ProcessTransitions.
  // Lexicographical order, with comparing the state only for "Element".

  class PairComparator {
   public:
    inline bool operator () (const pair<Label, Element> &p1, const pair<Label, Element> &p2) {
      if (p1.first < p2.first) return true;
      else if (p1.first > p2.first) return false;
      else {
        return p1.second.state < p2.second.state;
      }
    }
  };


  // ProcessTransitions handles transitions out of this subset of states.
  // Ignores epsilon transitions (epsilon closure already handled that).
  // Does not consider final states.  Breaks the transitions up by ilabel,
  // and creates a new transition in determinized FST, for each ilabel.
  // Does this by creating a big vector of pairs <Label, Element> and then sorting them
  // using a lexicographical ordering, and calling ProcessTransition for each range
  // with the same ilabel.
  // Side effects on repository, and (via ProcessTransition) on Q_, hash_,
  // and output_arcs_.
  void ProcessTransitions(const vector<Element> &closed_subset, OutputStateId state) {
    vector<pair<Label, Element> > all_elems;
    {  // Push back into "all_elems", elements corresponding to all non-epsilon-input transitions
      // out of all states in "closed_subset".
      typename vector<Element>::const_iterator iter = closed_subset.begin(),
          end = closed_subset.end();
      for (; iter != end; ++iter) {
        const Element &elem = *iter;
        for (ArcIterator<Fst<Arc> > aiter(*ifst_, elem.state);
             !aiter.Done(); aiter.Next()) {
          const Arc &arc = aiter.Value();
          if (arc.ilabel != 0) {  // Non-epsilon transition -- ignore epsilons here.
            pair<Label, Element> this_pr;
            this_pr.first = arc.ilabel;
            Element &next_elem(this_pr.second);
            next_elem.state = arc.nextstate;
            next_elem.weight = Times(elem.weight, arc.weight);
            if (arc.olabel == 0) // output epsilon-- this is simple case so
                                 // handle separately for efficiency
              next_elem.string = elem.string;
            else {
              vector<Label> seq;
              repository_.SeqOfId(elem.string, &seq);
              seq.push_back(arc.olabel);
              next_elem.string = repository_.IdOfSeq(seq);
            }
            all_elems.push_back(this_pr);
          }
        }
      }
    }
    PairComparator pc;
    std::sort(all_elems.begin(), all_elems.end(), pc);
    // now sorted first on input label, then on state.
    typedef typename vector<pair<Label, Element> >::const_iterator PairIter;
    PairIter cur = all_elems.begin(), end = all_elems.end();
    vector<Element> this_subset;
    while (cur != end) {
      // Process ranges that share the same input symbol.
      Label ilabel = cur->first;
      this_subset.clear();
      while (cur != end && cur->first == ilabel) {
        this_subset.push_back(cur->second);
        cur++;
      }
      // We now have a subset for this ilabel.
      ProcessTransition(state, ilabel, &this_subset);
    }
  }

  // SubsetToStateId converts a subset (vector of Elements) to a StateId in the output
  // fst.  This is a hash lookup; if no such state exists, it adds a new state to the hash
  // and adds a new pair to the queue.
  // Side effects on hash_ and Q_, and on output_arcs_ [just affects the size].
  OutputStateId SubsetToStateId(const vector<Element> &subset) {  // may add the subset to the queue.
    typedef typename SubsetHash::iterator IterType;
    IterType iter = hash_.find(&subset);
    if (iter == hash_.end()) {  // was not there.
      vector<Element> *new_subset = new vector<Element>(subset);
      OutputStateId new_state_id = (OutputStateId) output_arcs_.size();
      bool ans = hash_.insert(std::pair<const vector<Element>*,
                                        OutputStateId>(new_subset,
                                                       new_state_id)).second;
      assert(ans);
      output_arcs_.push_back(vector<TempArc>());
      if (allow_partial_ == false) {
        // If --allow-partial is not requested, we do the old way.
        Q_.push_front(pair<vector<Element>*, OutputStateId>(new_subset,  new_state_id));
      } else {
        // If --allow-partial is requested, we do breadth first search. This
        // ensures that when we return partial results, we return the states
        // that are reachable by the fewest steps from the start state.
        Q_.push_back(pair<vector<Element>*, OutputStateId>(new_subset,  new_state_id));
      }
      return new_state_id;
    } else {
      return iter->second;  // the OutputStateId.
    }
  }


  // ProcessSubset does the processing of a determinized state, i.e. it creates
  // transitions out of it and adds new determinized states to the queue if necessary.
  // The first stage is "EpsilonClosure" (follow epsilons to get a possibly larger set
  // of (states, weights)).  After that we ignore epsilons.  We process the final-weight
  // of the state, and then handle transitions out (this may add more determinized states
  // to the queue).
  void ProcessSubset(const pair<vector<Element>*, OutputStateId> & pair) {
    const vector<Element> *subset = pair.first;
    OutputStateId state = pair.second;

    vector<Element> closed_subset;  // subset after epsilon closure.
    epsilon_closure_.GetEpsilonClosure(*subset, &closed_subset);

    // Now follow non-epsilon arcs [and also process final states]
    ProcessFinal(closed_subset, state);

    // Now handle transitions out of these states.
    ProcessTransitions(closed_subset, state);
  }

  void Debug();

  DISALLOW_COPY_AND_ASSIGN(DeterminizerStar);
  deque<pair<vector<Element>*, OutputStateId> > Q_;  // queue of subsets to be processed.

  vector<vector<TempArc> > output_arcs_;  // essentially an FST in our format.

  const Fst<Arc> *ifst_;
  float delta_;
  int max_states_;
  bool determinized_; // used to check usage.
  bool allow_partial_;  // output paritial results or not
  bool is_partial_;     // if we get partial results or not
  SubsetKey hasher_;  // object that computes keys-- has no data members.
  SubsetEqual equal_;  // object that compares subsets-- only data member is delta_.
  SubsetHash hash_;  // hash from Subset to StateId in final Fst.

  StringRepository<Label, StringId> repository_;  // associate integer id's with sequences of labels.
  EpsilonClosure epsilon_closure_;
};


template<class F>
bool DeterminizeStar(F &ifst, MutableFst<typename F::Arc> *ofst,
                     float delta, bool *debug_ptr, int max_states,
                     bool allow_partial) {
  ofst->SetOutputSymbols(ifst.OutputSymbols());
  ofst->SetInputSymbols(ifst.InputSymbols());
  DeterminizerStar<F> det(ifst, delta, max_states, allow_partial);
  det.Determinize(debug_ptr);
  det.Output(ofst);
  return det.IsPartial();
}


template<class F>
bool DeterminizeStar(F &ifst,
                     MutableFst<GallicArc<typename F::Arc> > *ofst, float delta,
                     bool *debug_ptr, int max_states,
                     bool allow_partial) {
  ofst->SetOutputSymbols(ifst.InputSymbols());
  ofst->SetInputSymbols(ifst.InputSymbols());
  DeterminizerStar<F> det(ifst, delta, max_states, allow_partial);
  det.Determinize(debug_ptr);
  det.Output(ofst);
  return det.IsPartial();
}

template<class F>
void DeterminizerStar<F>::EpsilonClosure::
            GetEpsilonClosure(const vector<Element> &input_subset,
                                       vector<Element> *output_subset) {
  ecinfo_.resize(0);
  size_t size = input_subset.size();
  // find whether input fst is known to be sorted in input label.
  bool sorted =
          ((ifst_->Properties(kILabelSorted, false) & kILabelSorted) != 0);

  // size is still the input_subset.size()
  for (size_t i = 0; i < size; i++) {
    ExpandOneElement(input_subset[i], sorted, input_subset[i].weight, true);
  }

  size_t s = queue_2_.size();
  if (s == 0) {
    *output_subset = input_subset;
    return;
  } else {
    // queue_2 not empty. Need to create the vector<info>
    for (size_t i = 0; i < size; i++) {
      // the weight has not been processed yet,
      // so put all of them in the "weight_to_process"
      ecinfo_.push_back(EpsilonClosureInfo(input_subset[i],
                                           input_subset[i].weight,
                                           false));
      ecinfo_.back().element.weight = Weight::Zero(); // clear the weight

      if (id_to_index_.size() < input_subset[i].state + 1) {
        id_to_index_.resize(2 * input_subset[i].state + 1, -1);
      }
      id_to_index_[input_subset[i].state] = ecinfo_.size() - 1;
    }
  }

  {
    Element elem;
    elem.weight = Weight::Zero();
    for (size_t i = 0; i < s; i++) {
      elem.state = queue_2_[i].state;
      elem.string = queue_2_[i].string;
      AddOneElement(elem, queue_2_[i].weight);
    }
    queue_2_.resize(0);
  }

  int counter = 0; // relates to max-states option, used for test.
  while (!queue_.empty()) {
    InputStateId id = queue_.front();

    // no need to check validity of the index
    // since anything in the queue we are sure they're in the "virtual set"
    int index = id_to_index_[id];
    EpsilonClosureInfo &info = ecinfo_[index];
    Element &elem = info.element;
    Weight unprocessed_weight = info.weight_to_process;

    elem.weight = Plus(elem.weight, unprocessed_weight);
    info.weight_to_process = Weight::Zero();

    info.in_queue = false;
    queue_.pop_front();

    if (max_states_ > 0 && counter++ > max_states_) {
      std::cerr << "Determinization aborted since looped more than "
                << max_states_ << " times during epsilon closure.\n";
      throw std::runtime_error("looped more than max-states times"
         " in determinization");
    }

    // generally we need to be careful about iterator-invalidation problem
    // here we pass a reference (elem), which could be an issue.
    // In the beginning of ExpandOneElement, we make a copy of elem.string
    // to avoid that issue
    ExpandOneElement(elem, sorted, unprocessed_weight);
  }

  {
    // this sorting is based on StateId
    sort(ecinfo_.begin(), ecinfo_.end());

    output_subset->clear();

    size = ecinfo_.size();
    output_subset->reserve(size);
    for (size_t i = 0; i < size; i++) {
      EpsilonClosureInfo& info = ecinfo_[i];
      if (info.weight_to_process != Weight::Zero()) {
        info.element.weight = Plus(info.element.weight, info.weight_to_process);
      }
      output_subset->push_back(info.element);
    }
  }
}

template<class F>
void DeterminizerStar<F>::EpsilonClosure::
     AddOneElement(const Element &elem, const Weight &unprocessed_weight) {
  // first we try to find the element info in the ecinfo_ vector
  int index = -1;
  if (elem.state < id_to_index_.size()) {
    index = id_to_index_[elem.state];
  }
  if (index != -1) {
    if (index >= ecinfo_.size()) {
      index = -1;
    }
    // since ecinfo_ might store outdated information, we need to check
    else if (ecinfo_[index].element.state != elem.state) {
      index = -1;
    }
  }

  if (index == -1) {
    // was no such StateId: insert and add to queue.
    ecinfo_.push_back(EpsilonClosureInfo(elem, unprocessed_weight, true));
    size_t size = id_to_index_.size();
    if (size < elem.state + 1) {
      // double the size to reduce memory operations
      id_to_index_.resize(2 * elem.state + 1, -1);
    }
    id_to_index_[elem.state] = ecinfo_.size() - 1;
    queue_.push_back(elem.state);

  } else {  // one is already there.  Add weights.
    EpsilonClosureInfo &info = ecinfo_[index];
    if (info.element.string != elem.string) {
      std::cerr << "DeterminizerStar: FST was not functional "
        "-> not determinizable\n";
      { // Print some debugging information.  Can be helpful to debug
        // the inputs when FSTs are mysteriously non-functional.
        vector<Label> tmp_seq;
        repository_->SeqOfId(info.element.string, &tmp_seq);
        std::cerr << "First string: ";
        for (size_t i = 0; i < tmp_seq.size(); i++)
          std::cerr << tmp_seq[i] << " ";
        std::cerr << "\nSecond string: ";
        repository_->SeqOfId(elem.string, &tmp_seq);
        for (size_t i = 0; i < tmp_seq.size(); i++)
          std::cerr << tmp_seq[i] << " ";
        std::cerr << "\n";
      }
      throw std::runtime_error("Non-functional FST: cannot determinize.\n");
    }

    info.weight_to_process =
          Plus(info.weight_to_process, unprocessed_weight);

    if (!info.in_queue) {
      // this is because the code in "else" below: the 
      // iter->second.weight_to_process might not be Zero()
      Weight weight = Plus(info.element.weight, info.weight_to_process);

      // What is done below is, we propagate the weight (by adding them
      // to the queue only when the change is big enough;
      // otherwise we just store the weight, until before returning
      // we add the element.weight and weight_to_process together
      if (! ApproxEqual(weight, info.element.weight, delta_)) {
        // add extra part of weight to queue.
        info.in_queue = true;
        queue_.push_back(elem.state);
      }
    }
  }
}

template<class F>
void DeterminizerStar<F>::EpsilonClosure::ExpandOneElement(
                                          const Element &elem,
                                          bool sorted,
                                          const Weight &unprocessed_weight,
                                          bool save_to_queue_2) {
  StringId str = elem.string; // copy it here because there is an iterator-
                // - invalidation problem (it really happens for some FSTs)

  // now we are going to propagate the "unprocessed_weight"
  for (ArcIterator<Fst<Arc> > aiter(*ifst_, elem.state);
       !aiter.Done(); aiter.Next()) {
    const Arc &arc = aiter.Value();
    if (sorted && arc.ilabel > 0) {
      break;
      // Break from the loop: due to sorting there will be no
      // more transitions with epsilons as input labels.
  
