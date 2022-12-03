
/**
 * THIS IS SCHOOL PROJECT FOR SUBJECT: AUTOMATA AND GRAMMARS
 * AUTHOR: KRYSTOF KLEN
 */

#ifndef __PROGTEST__
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <variant>
#include <vector>

using State = unsigned int;
using Symbol = uint8_t;

struct NFA {
    std::set<State> m_States;
    std::set<Symbol> m_Alphabet;
    std::map<std::pair<State, Symbol>, std::set<State>> m_Transitions;
    State m_InitialState;
    std::set<State> m_FinalStates;
};

struct DFA {
    std::set<State> m_States;
    std::set<Symbol> m_Alphabet;
    std::map<std::pair<State, Symbol>, State> m_Transitions;
    State m_InitialState;
    std::set<State> m_FinalStates;
};
#endif
//______________________________________________________________________________________________________________________
const size_t DEAD_NODE = 100000; // No automata state can have key == DEAD_NODE!!!
void print(const NFA & nfa, const std::string & msg){
    std::cout<<"****************    "<<msg<<"    ********************\n";
    for(State st : nfa.m_States){
        std::cout<<st<<" | ";
        for(Symbol  sy : nfa.m_Alphabet){
            std::cout<<"          "<<sy<<" -> { ";
            auto it = nfa.m_Transitions.find({st,sy});
            if(it != nfa.m_Transitions.end()){
                for(auto target : it->second){
                    std::cout<<target<<" ";
                }
            }
            std::cout<<" }";
        }
        std::cout<<"\n";
    }
    std::cout<<"Itinial state: "<<nfa.m_InitialState<<std::endl;
    std::cout<<"Final states: ";
    for(auto x : nfa.m_FinalStates) std::cout<<x<<" ";
    std::cout<<"\n";
}
void print(DFA & dfa, std::deque<State> group, const std::string & msg){
    std::cout<<"****************    "<<msg<<"    ********************\n";

    for(State st : dfa.m_States){
        State gs = group[st];
        std::cout<<gs<<" | ";
        for(Symbol  sy : dfa.m_Alphabet){
            std::cout<<"          "<<sy<<" -> { ";
            auto it = dfa.m_Transitions.find({st,sy});
            if(it != dfa.m_Transitions.end()){
                if(it->second == DEAD_NODE) {
                    std::cout<<DEAD_NODE<<" ";
                }else{
                    State groupState = group[it->second];
                    std::cout<<groupState<<" ";
                }
            }
            std::cout<<" }";
        }
        std::cout<<"\n";
    }
    std::cout<<"Itinial state: "<<dfa.m_InitialState<<std::endl;
    std::cout<<"Final states: ";
    for(auto x : dfa.m_FinalStates) std::cout<<x<<" ";
    std::cout<<"\n";
}
void printS(NFA & nfa, const std::string & msg){
    std::cout<<"NFA ";
    for(Symbol sy : nfa.m_Alphabet){
        std::cout<<sy<<" ";
    }
    std::cout<<"\n";
    for(State st : nfa.m_States){
        if(nfa.m_FinalStates.find(st) != nfa.m_FinalStates.end()){
            std::cout<<"<";
        }
        if(nfa.m_InitialState == st ){
            std::cout<<">";
        }
        std::cout<<" "<<st;
        for(Symbol sy : nfa.m_Alphabet){
            auto tr = nfa.m_Transitions.find({st,sy});
            if(tr == nfa.m_Transitions.end()) {
                std::cout<<" -";
            }else{
                std::cout<<" ";
                size_t x = 0;
                for(auto to : tr->second){
                    std::cout<<to;
                    x++;
                    if(x<tr->second.size()){
                        std::cout<<"|";
                    }
                }
            }
        }
        std::cout<<"\n";
    }

}
void printS(DFA & dfa, const std::string & msg){
    std::cout<<"________________ "<< msg<<std::endl;
    std::cout<<"DFA ";
    for(Symbol sy : dfa.m_Alphabet){
        std::cout<<sy<<" ";
    }
    std::cout<<"\n";
    for(State st : dfa.m_States){
        if(dfa.m_FinalStates.find(st) != dfa.m_FinalStates.end()){
            std::cout<<"<";
        }
        if(dfa.m_InitialState == st ){
            std::cout<<">";
        }
        std::cout<<" "<<st;
        for(Symbol sy : dfa.m_Alphabet){
            auto tr = dfa.m_Transitions.find({st,sy});
            if(tr == dfa.m_Transitions.end()) {
                std::cout<<" -";
            }else{
                std::cout<< " "<<tr->second;
            }
        }
        std::cout<<"\n";
    }

}

void normalize(DFA & dfa){
    std::map<State,State> dct;
    std::set<State> final;
    for(State st : dfa.m_States){
        State cr = dct.size();
        dct.insert({st,cr});
        if(dfa.m_FinalStates.find(st) != dfa.m_FinalStates.end()){
            final.insert(cr);
        }
    }
    std::map<std::pair<State, Symbol>, State> tr;
    for(auto trans : dfa.m_Transitions){
        State from = dct[trans.first.first];
        State to = dct[trans.second];
        Symbol sy = trans.first.second;
        tr.insert({{from,sy},to});
    }
    dfa.m_Transitions.clear();
    for(auto trans :tr){
        dfa.m_Transitions.insert(trans);
    }
    dfa.m_States.clear();
    for(auto st : dct){
        dfa.m_States.insert(st.second);
    }
    dfa.m_FinalStates.clear();
    dfa.m_FinalStates = final;

}
//______________________________________________________________________________________________________________________
/**
 * Constructs DFA for empty language as teacher's instructions
 * @param res
 * @param alphabet
 */
void constructAutotamaGeneratingEmptyL( DFA & res, const std::set<Symbol> & alphabet ){
    res.m_States.clear();
    res.m_Transitions.clear();
    res.m_FinalStates.clear();

    res.m_States.insert(0);
    res.m_InitialState = 0;
    for(Symbol sy : alphabet){
        res.m_Alphabet.insert(sy);
    }
}
//______________________________________________________________________________________________________________________
void copyAlbhabet(const std::set<Symbol> & from, std::set<Symbol> & to){
    for(Symbol sy : from){
        to.insert(sy);
    }
}
//______________________________________________________________________________________________________________________
/**
 * Helper method for mergeAutomata.
 * Adds states autoata can reach from starting point
 * after removing epsilon transitions from initial state
 * @param a
 * @param result
 * @param mappingResult = state from automata a mapped to state in the merged automata
 */
void addTransitionFromInitState(const NFA & a, NFA & result, std::map<State,State> & mappingResult ){
    // add transitions from new initial state to states where I can get from original initial state
    for(Symbol sy : a.m_Alphabet){
        auto it = a.m_Transitions.find({a.m_InitialState,sy});
        if( it != a.m_Transitions.end()){
            // there are transitions from the original init state
            for(State stReacheble : it->second){
                // add transition from init to reachable state
                State stInResult = mappingResult.at(stReacheble);
                result.m_Transitions[{result.m_InitialState,sy}].insert(stInResult);
            }
        }
    }
}
//______________________________________________________________________________________________________________________
void copyTransitions(const NFA & from ,NFA & to, std::map<State,State> & aResult){
    for(auto tr : from.m_Transitions){
        State stFrom = aResult.at(tr.first.first);
        Symbol sy = tr.first.second;
        auto it = to.m_Transitions.find({ stFrom, sy });
        if(it == to.m_Transitions.end()){
            to.m_Transitions.insert({{stFrom,sy},{}});
        }
        for(State x : tr.second){
            State stTo = aResult.at(x);
            to.m_Transitions.at({stFrom,sy}).insert(stTo);
        }
    }
}
//______________________________________________________________________________________________________________________
void setFinalStates(const NFA & a, NFA & result, std::map<State,State> & aResult){
    for(State st : a.m_FinalStates){
        State stRes = aResult.at(st);
        result.m_FinalStates.insert(stRes);
        if(st == a.m_InitialState){
            result.m_FinalStates.insert(result.m_InitialState);
        }
    }
}
//______________________________________________________________________________________________________________________
/**
 * Makes unifications of 2 Nondeterministic finite state machines.
 * Creates a new initial state with with transitions to states
 * withins it's epsilon closure.
 * @param a
 * @param b
 * @param result
 */
void mergeAutomata(const NFA & a, const NFA & b, NFA & result){
    std::set<State> reachableFromStartA, reachableFromStartB;
    std::map<State,State> aResult, bResult;

    // copy all states with
    for(State st : a.m_States){
        State cState = result.m_States.size() + 1;
        aResult.insert({st,cState});
        result.m_States.insert(cState);
    }
    for(State st : b.m_States){
        State cState = result.m_States.size() + 1;
        bResult.insert({st,cState});
        result.m_States.insert(cState);
    }

    State initState = 0;
    result.m_States.insert(initState);
    result.m_InitialState = initState;

    // add alphabet
    copyAlbhabet(a.m_Alphabet,result.m_Alphabet);
    copyAlbhabet(b.m_Alphabet,result.m_Alphabet);

    addTransitionFromInitState(a,result,aResult);
    addTransitionFromInitState(b,result,bResult);

    // copy all transitions
    copyTransitions(a,result,aResult);
    copyTransitions(b,result,bResult);

    // set final states
    setFinalStates(a,result,aResult);
    setFinalStates(b,result,bResult);
}
//______________________________________________________________________________________________________________________
void determinize(const NFA& nfa, DFA& dfa){
    copyAlbhabet(nfa.m_Alphabet,dfa.m_Alphabet);
   // dfa.m_InitialState = 0;
    using set_state = std::set<State>;
    std::map<set_state,State> mapping;
    std::queue<set_state> q;

    q.push({nfa.m_InitialState});
    mapping.insert({{nfa.m_InitialState},mapping.size()});

    while(!q.empty()){
        auto curr = q.front();
        State from = mapping.at(curr);
        dfa.m_States.insert(from);
        for(Symbol sy : dfa.m_Alphabet){
            set_state next;
            for(State sFrom : curr){
                if(nfa.m_FinalStates.find(sFrom) != nfa.m_FinalStates.end()){
                    dfa.m_FinalStates.insert(from);
                }
                auto trans = nfa.m_Transitions.find({sFrom,sy});
                if(trans == nfa.m_Transitions.end()) continue;

                for(auto sTo : trans->second){
                    next.insert(sTo);
                }
            }
            if( !next.empty() && mapping.find(next) == mapping.end() ){
                q.push(next);
                mapping.insert({next,mapping.size()});
            }
            if(!next.empty()){
                State to = mapping.find(next)->second;
                dfa.m_Transitions.insert({{from,sy},to});
            }
        }

        q.pop();
    }
    dfa.m_InitialState = mapping[{nfa.m_InitialState}];
}
//______________________________________________________________________________________________________________________
void reverseTransitions(const std::map<std::pair<State, Symbol>, State> & transitions,
                        std::map<std::pair<State, Symbol>, std::set<State>> & converted ){
    for(auto trans : transitions){
        State from  = trans.second;
        State to = trans.first.first;
        Symbol sy = trans.first.second;
        converted[{from,sy}].insert(to);
    }
}
//______________________________________________________________________________________________________________________
bool removeUseless(DFA & dfa){
    std::map<std::pair<State, Symbol>, std::set<State>>  converted;
    reverseTransitions(dfa.m_Transitions, converted);
    std::set<State> useful;

    for(State fSt : dfa.m_FinalStates){
        // final states useful by definition
        useful.insert(fSt);
    }
    for(State fSt : dfa.m_FinalStates){
        // flooding algorithm
        std::queue<State> q;
        q.push(fSt);
        while( ! q.empty()){
            auto curr = q.front();
            for(Symbol sy : dfa.m_Alphabet){
                auto trans = converted.find({curr,sy});
                if(trans == converted.end()){
                    continue;
                }
                for( auto to : trans->second){
                    if( useful.find(to) != useful.end() ){
                        // State already closed
                        continue;
                    }
                    q.push(to);
                    useful.insert(to);
                }
            }
            q.pop();
        }
    }

    // load useful states to DFA, final states remain the same
    dfa.m_States.clear();
    for(State uSt : useful){
        dfa.m_States.insert(uSt);
    }

    if(dfa.m_States.find(dfa.m_InitialState) == dfa.m_States.end()){
        // automata accepts empty language
        constructAutotamaGeneratingEmptyL(dfa,dfa.m_Alphabet);
        return false;
    }

    // update transitions
    for(auto it = dfa.m_Transitions.begin(); it != dfa.m_Transitions.end();){
        State from = it->first.first;
        State to = it->second;

        if(dfa.m_States.find(from) == dfa.m_States.end()
        || dfa.m_States.find(to) == dfa.m_States.end()){
            dfa.m_Transitions.erase(it++);
        }else{
            it++;
        }
    }
    return true;
}
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
void cartesianProductStates(const std::set<State> &a,
                            const std::set<State> &b,
                            std::map < std::pair<State,State>, State > & carP){
    for(State as : a){
        for(State bs : b){
            carP.insert({{as,bs},carP.size()});
        }
    }
}
//______________________________________________________________________________________________________________________

bool finalAtIntersection(State a, State b,
                         const std::set<State> & finalStatesA,
                         const std::set<State> & finalStatesB ){
    return
        finalStatesA.find(a) != finalStatesA.end()
        &&
        finalStatesB.find(b) != finalStatesB.end();
}
//______________________________________________________________________________________________________________________
/**
 * Makes an intersection of 2 Nondeterministic finite state machines
 * using cartesian product approach.
 * @param a
 * @param b
 * @param result
 */
void automataIntersectionCartesianProduct(const NFA& a, const NFA& b, NFA & result){
    std::map < std::pair<State,State>, State > carP;
    cartesianProductStates(a.m_States,b.m_States,carP);
    copyAlbhabet(a.m_Alphabet,result.m_Alphabet);
    copyAlbhabet(b.m_Alphabet,result.m_Alphabet);


    result.m_InitialState = carP[{a.m_InitialState,b.m_InitialState}];

    for(auto carProdSt : carP){
        State fA = carProdSt.first.first;
        State fB = carProdSt.first.second;
        result.m_States.insert(carProdSt.second);
        if(finalAtIntersection(fA,fB,a.m_FinalStates,b.m_FinalStates) ){
            result.m_FinalStates.insert(carProdSt.second);
        }

        for(Symbol sy : result.m_Alphabet){
            bool aOK = false;
            bool bOK = false;
            auto transA = a.m_Transitions.find({fA,sy});
            if(transA != a.m_Transitions.end()){
                aOK = true;
            }
            auto transB = b.m_Transitions.find({fB,sy});
            if(transB != b.m_Transitions.end()){
                bOK = true;
            }

            if(!aOK || !bOK) continue;

            for(auto toA  : transA->second){
                for(auto toB  : transB->second){
                    State to = carP[{toA,toB}];
                    result.m_Transitions[{carProdSt.second,sy}].insert(to);
                }
            }


        }
    }

}
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
void splitIntoFinishNonFinish(DFA & orig, std::deque<State> & group){
    State iNonFinish = 0;
    State iFinish = 0;
    bool iNonFinishSet = false;
    bool iFinishSet = false;

    for(State st = 0; st<group.size();st++){
        if(orig.m_FinalStates.find(st) != orig.m_FinalStates.end()){

            if(!iFinishSet){
                iFinish = st;
                iFinishSet = true;
            }
            group[st] = iFinish;
        } else {

            if( !iNonFinishSet){
                iNonFinish = st;
                iNonFinishSet = true;
            }
            group[st] = iNonFinish;
        }
    }
}
//______________________________________________________________________________________________________________________
/**
 * Makes next partitioning in the process of minimalizing DFA
 * @param orig
 * @param group = mapping states -> groups
 * @return
 */
bool nextStep(const DFA & orig, std::deque<State> & group){
    bool changeOccured = false;
    std::deque<State> prevGroud = group;

    // allocate new table for better orientation
    State** curr = new State*[group.size()];
    for(State i = 0; i<group.size(); i++){
        curr[i] = new State[orig.m_Alphabet.size()];
    }

    // fill table
    for(State st = 0; st<group.size(); st++){
        State syIndx = 0;
        for(Symbol sy : orig.m_Alphabet){
            auto tr = orig.m_Transitions.find({st,sy});
            if(tr == orig.m_Transitions.end()) {
                curr[st][syIndx] = DEAD_NODE;
                syIndx++;
                continue;
            }
            curr[st][syIndx] = group[tr->second];
            syIndx++;
        }
    }

    std::map<State,std::vector<State>> nextRound;
    // first iteration -> insert groups with line
    for(State i = 0; i<group.size();i++) {
        // iterate states
        State grpFrom = group[i];
        std::vector<State> vcLine;
        for (State j = 0; j < orig.m_Alphabet.size(); j++) {
           // if (curr[i][j] == DEAD_NODE) continue;
            vcLine.push_back(curr[i][j]);
        }
        // insert current groups
        auto it = nextRound.find(grpFrom);
        if (it == nextRound.end() ) {
            nextRound.insert({grpFrom, vcLine});
        }
    }
    // second iteration -> compare lines
    std::set<State> newGroups;
    for(State i = 0; i<group.size(); i++){
        State grpFrom = group[i];
        const std::vector<State> & antcLine = nextRound.at(grpFrom);
        std::vector<State> currLine;
        // get curr line
        for (State j = 0; j < orig.m_Alphabet.size(); j++){
            //if(curr[i][j]==DEAD_NODE) continue;
            currLine.push_back(curr[i][j]);
        }

        // line does not match
        if( currLine != antcLine ){
            changeOccured = true;
            // check new groups
            bool matchFound = false;
            for(State ng : newGroups){
                auto found = nextRound.find(ng);
                const std::vector<State> & possible = found->second;
                if( prevGroud[found->first]==group[i] && currLine == possible ){
                    matchFound = true;
                    group[i] = ng;
                }
            }
            if( !matchFound){
                // completely new group
                nextRound.insert({i,currLine});
                newGroups.insert(i);
                group[i] = i;
            }
        }
    }
    // deallocate memory
    for(State i = 0; i<group.size(); i++){
        delete []  curr[i];
    }
    delete [] curr;
    return changeOccured;
}
//______________________________________________________________________________________________________________________
/**
 * Assembles a new DFA after minimalization process
 * @param orig
 * @param group
 * @param output
 */
void createAutomata(DFA & orig, std::deque<State> & group, DFA & output ){

    output.m_InitialState = group[orig.m_InitialState];

    copyAlbhabet(orig.m_Alphabet,output.m_Alphabet);

   // copy states
    for(State st  :orig.m_States){
        output.m_States.insert(group[st]);
    }

    for(State st : output.m_States){
        for(Symbol sy  : orig.m_Alphabet){
            auto trans = orig.m_Transitions.find({st,sy});
            if(trans != orig.m_Transitions.end()){
                output.m_Transitions.insert({ { trans->first.first,sy   }   , group[trans->second]});
            }
        }
    }
    for(State st  : orig.m_FinalStates){
        output.m_FinalStates.insert(group[st]);
    }
    normalize(output);
}
//______________________________________________________________________________________________________________________
/**
 * Finds the smallest possible determistic finite state machine
 * for language accepted by @param dfa
 * @param dfa
 * @param resultFINAL
 */
void minimalize(DFA & dfa, DFA & resultFINAL){
    std::deque<State> groups(dfa.m_States.size());
    bool acceptsNONemptyL = removeUseless(dfa);
    if(!acceptsNONemptyL){
        constructAutotamaGeneratingEmptyL(resultFINAL,dfa.m_Alphabet);
        return;
    }
    normalize(dfa);
    splitIntoFinishNonFinish(dfa,groups);
    while(nextStep(dfa,groups)){
    }
    createAutomata(dfa,groups,resultFINAL);
}
//______________________________________________________________________________________________________________________


//______________________________________________________________________________________________________________________
/**
 * Makes  intersection of 2 NFAs
 * @param a
 * @param b
 * @return
 */
DFA intersect(const NFA& a, const NFA& b){
    DFA resultDFA;

    NFA resultNFA;
    automataIntersectionCartesianProduct(a,b,resultNFA);
    DFA resDFA;
    determinize(resultNFA,resDFA);
    minimalize(resDFA,resultDFA);
    return resultDFA;
}
/**
 * Makes unification of 2 NFAs
 * @param a
 * @param b
 * @return
 */
DFA unify(const NFA& a, const NFA& b) {
    DFA resultDFA;

    NFA resultNFA;
    mergeAutomata(a, b, resultNFA);
    DFA resDFA;
    determinize(resultNFA, resDFA);
    minimalize(resDFA,resultDFA);
    return resultDFA;
}


#ifndef __PROGTEST__
//________________________________________________________________________________________________
bool operator==(const DFA& a, const DFA& b)
{
    return std::tie(a.m_States, a.m_Alphabet, a.m_Transitions, a.m_InitialState, a.m_FinalStates) == std::tie(b.m_States, b.m_Alphabet, b.m_Transitions, b.m_InitialState, b.m_FinalStates);
}
void testIntersect(){
    NFA a1{
            {0, 1, 2},
            {'a', 'b'},
            {
                    {{0, 'a'}, {0, 1}},
                    {{0, 'b'}, {0}},
                    {{1, 'a'}, {2}},
            },
            0,
            {2},
    };
    NFA a2{
            {0, 1, 2},
            {'a', 'b'},
            {
                    {{0, 'a'}, {1}},
                    {{1, 'a'}, {2}},
                    {{2, 'a'}, {2}},
                    {{2, 'b'}, {2}},
            },
            0,
            {2},
    };
    DFA a{
            {0, 1, 2, 3, 4},
            {'a', 'b'},
            {
                    {{0, 'a'}, {1}},
                    {{1, 'a'}, {2}},
                    {{2, 'a'}, {2}},
                    {{2, 'b'}, {3}},
                    {{3, 'a'}, {4}},
                    {{3, 'b'}, {3}},
                    {{4, 'a'}, {2}},
                    {{4, 'b'}, {3}},
            },
            0,
            {2},
    };
    assert(intersect(a1, a2) == a);
    NFA d1{
            {0, 1, 2, 3},
            {'i', 'k', 'q'},
            {
                    {{0, 'i'}, {2}},
                    {{0, 'k'}, {1, 2, 3}},
                    {{0, 'q'}, {0, 3}},
                    {{1, 'i'}, {1}},
                    {{1, 'k'}, {0}},
                    {{1, 'q'}, {1, 2, 3}},
                    {{2, 'i'}, {0, 2}},
                    {{3, 'i'}, {3}},
                    {{3, 'k'}, {1, 2}},
            },
            0,
            {2, 3},
    };
    NFA d2{
            {0, 1, 2, 3},
            {'i', 'k'},
            {
                    {{0, 'i'}, {3}},
                    {{0, 'k'}, {1, 2, 3}},
                    {{1, 'k'}, {2}},
                    {{2, 'i'}, {0, 1, 3}},
                    {{2, 'k'}, {0, 1}},
            },
            0,
            {2, 3},
    };
    DFA d{
            {0, 1, 2, 3},
            {'i', 'k', 'q'},
            {
                    {{0, 'i'}, {1}},
                    {{0, 'k'}, {2}},
                    {{2, 'i'}, {3}},
                    {{2, 'k'}, {2}},
                    {{3, 'i'}, {1}},
                    {{3, 'k'}, {2}},
            },
            0,
            {1, 2, 3},
    };

    DFA r = intersect(d1,d2);
    assert(r==d);
    NFA c1{
            {0, 1, 2, 3, 4},
            {'a', 'b'},
            {
                    {{0, 'a'}, {1}},
                    {{0, 'b'}, {2}},
                    {{2, 'a'}, {2, 3}},
                    {{2, 'b'}, {2}},
                    {{3, 'a'}, {4}},
            },
            0,
            {1, 4},
    };
    NFA c2{
            {0, 1, 2},
            {'a', 'b'},
            {
                    {{0, 'a'}, {0}},
                    {{0, 'b'}, {0, 1}},
                    {{1, 'b'}, {2}},
            },
            0,
            {2},
    };
    DFA c{
            {0},
            {'a', 'b'},
            {},
            0,
            {},
    };
    assert(intersect(c1, c2) == c);
}
void testUnion(){
    NFA b1{
            {0, 1, 2, 3, 4},
            {'a', 'b'},
            {
                    {{0, 'a'}, {1}},
                    {{0, 'b'}, {2}},
                    {{2, 'a'}, {2, 3}},
                    {{2, 'b'}, {2}},
                    {{3, 'a'}, {4}},
            },
            0,
            {1, 4},
    };
    NFA b2{
            {0, 1, 2, 3, 4},
            {'a', 'b'},
            {
                    {{0, 'b'}, {1}},
                    {{1, 'a'}, {2}},
                    {{2, 'b'}, {3}},
                    {{3, 'a'}, {4}},
                    {{4, 'a'}, {4}},
                    {{4, 'b'}, {4}},
            },
            0,
            {4},
    };
    DFA b{
            {0, 1, 2, 3, 4, 5, 6, 7, 8},
            {'a', 'b'},
            {
                    {{0, 'a'}, {1}},
                    {{0, 'b'}, {2}},
                    {{2, 'a'}, {3}},
                    {{2, 'b'}, {4}},
                    {{3, 'a'}, {5}},
                    {{3, 'b'}, {6}},
                    {{4, 'a'}, {7}},
                    {{4, 'b'}, {4}},
                    {{5, 'a'}, {5}},
                    {{5, 'b'}, {4}},
                    {{6, 'a'}, {8}},
                    {{6, 'b'}, {4}},
                    {{7, 'a'}, {5}},
                    {{7, 'b'}, {4}},
                    {{8, 'a'}, {8}},
                    {{8, 'b'}, {8}},
            },
            0,
            {1, 5, 8},
    };
    assert(unify(b1, b2) == b);
}

int main() {
    
    return 0;
}
#endif
