#include "Eigen/Dense"
#include "Msg.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

Tree::Tree(const Tree& t) {

    clone(t);
}

Tree::Tree(std::string newick){
    numTaxa = 0;
    std::vector<std::string> newickTokens = parseNewickString(newick);
    Node* p = nullptr;
    bool readingBl = false;
    for(int i = 0; i < newickTokens.size(); i++){
        std::string token = newickTokens[i];
        if(token == "("){
            Node* newNode = addNode();
            if(p == nullptr){
                root = newNode;
            }else{
                p->addNeighbor(newNode);
                newNode->addNeighbor(p);
                newNode->setAncestor(p);
            }
            p = newNode;
        }else if (token == ")" || token == ","){
            if(p->getAncestor() == nullptr)
                Msg::error("no anc found for p");
            p = p->getAncestor();
        }else if (token == ";"){
            if(p != root)
                Msg::error("expecting to be at root");
        }else if (token == ":"){
            readingBl = true;
        }else{
            if(readingBl == false){
                Node* newNode = addNode();
                newNode->addNeighbor(p);
                newNode->setAncestor(p);
                p->addNeighbor(newNode);
                newNode->setName(token);
                newNode->setIsTip(true);
                tips.push_back(newNode);
                numTaxa++;
                p = newNode;
            }else{
                double x = stod(token);
                if(p->getAncestor() != nullptr)
                    setBranch(p, p->getAncestor(), x);
                readingBl = false;
            }
        }
    }
    initializeDownPassSequence();
    int idx = numTaxa;
    int tIdx = 0;
    for (Node* p : downPassSequence)
        {
            if(p->getIsTip() == true)
                p->setIndex(tIdx++);
            if (p->getIsTip() == false)
                p->setIndex(idx++);
        }
}

Tree::~Tree(void) {

    deleteNodes();
}

Tree& Tree::operator=(const Tree& t) {

    if (this != &t)
        clone(t);
    return *this;
}

Node* Tree::addNode(void) {

    Node* newNode = new Node;
    newNode->setOffset((int)nodes.size());
    nodes.push_back(newNode);
    return newNode;
}

void Tree::checkBranchLengthsNeg(void){
    for(auto const& x : branchLengths)
        if(x.second < 0)
            Msg::error("negative branch lengths");
}

void Tree::clone(const Tree& t) {
    
    if (this->nodes.size() != t.nodes.size())
        {
        deleteNodes();
        for (int i=0; i<t.nodes.size(); i++)
            addNode();
        }
        
    this->numTaxa = t.numTaxa;
    this->root = this->nodes[t.root->getOffset()];
    
    for (int i=0; i<t.nodes.size(); i++)
        {
        Node* q = t.nodes[i];
        Node* p = this->nodes[i];
        p->setIndex(q->getIndex());
        p->setIsTip(q->getIsTip());
        p->setName(q->getName());
        if (q->getAncestor() != nullptr)
            p->setAncestor( this->nodes[q->getAncestor()->getOffset()] );
        else
            p->setAncestor(nullptr);
        std::set<Node*>& qNeighbors = q->getNeighbors();
        p->removeAllNeighbors();
        for (Node* qn : qNeighbors)
            p->addNeighbor( this->nodes[qn->getOffset()] );
        }
        
    initializeDownPassSequence();
    this->branchLengths.clear();
    for (BranchLengths::const_iterator it = t.branchLengths.begin(); it != t.branchLengths.end(); it++)
        {
        if(it->first.first == nullptr)
            Msg::error("stop");
        Node* e1 = this->nodes[it->first.first->getOffset()];
        Node* e2 = this->nodes[it->first.second->getOffset()];
        setBranch(e1, e2, it->second);
        }
}

void Tree::deleteNodes(void) {

    for (int i=0; i<nodes.size(); i++)
        delete nodes[i];
    nodes.clear();
}
    
double Tree::getBranchLength(Node* e1, Node* e2) {

    std::pair<Node*,Node*> key;
    initializeBranchLengthKey(key, e1, e2);
    BranchLengths::iterator it = branchLengths.find(key);
    if (it == branchLengths.end()){
        print("Tree with issue");
        std::cout << "for nodes: " << e1->getIndex() << " " << e2->getIndex() << std::endl;
        std::cout << "for nodes: " << e1 << " " << e2 << std::endl;
        std::cout << "e1 anc: " << e1->getAncestor()->getIndex() << std::endl;
        Msg::error("Could not find branch to return");
    }
    return it->second;
}

std::string Tree::getNewickString(void) {

    std::stringstream strm;
    writeTree(root, strm);
    strm << ";";
    return strm.str();
}

void Tree::initializeBranchLengthKey(std::pair<Node*,Node*>& key, Node* e1, Node* e2) {

    if (e1 < e2)
        key = std::make_pair(e1,e2);
    else
        key = std::make_pair(e2,e1);
}

void Tree::initializeDownPassSequence(void) {
    if(root == nullptr)
        Msg::error("root is nullptr");
    downPassSequence.clear();
    passDown(root, root);
    checkBranchLengthsNeg();
}

std::vector<std::string>  Tree::parseNewickString(std::string newickStr){
    
    std::vector<std::string> tokens;
    
    std::string token = "";
    for(int i = 0; i < newickStr.length(); i++){
        char c = newickStr[i];
        if(c == '(' || c == ')' || c==',' || c==':' ||c == ';'){
            if(token != ""){
                tokens.push_back(token);
                token = "";
            }
            tokens.push_back(std::string(1,c));
        }else{
            token += std::string(1,c);
        }
    }
    
    if(token != "")
        tokens.push_back(token);
    
//    for(int i = 0; i < tokens.size(); i++){
//        std::cout << i << " " << tokens[i] << std::endl;
//    }
//        
    
    return tokens;
}

void Tree::passDown(Node* p, Node* from) {

    if (p != nullptr)
        {
        std::set<Node*>& pNeighbors = p->getNeighbors();
        for (Node* d : pNeighbors)
            {
            if (d != from)
                passDown(d, p);
            }
        p->setAncestor(from);
        downPassSequence.push_back(p);
        }
}

void Tree::print(void) {
    
    showNode(root, 0);
    std::cout << "Postorder sequence: ";
    for (int i=0; i<downPassSequence.size(); i++)
        std::cout << downPassSequence[i]->getIndex() << " ";
    std::cout << std::endl;
    //printVCV();
}

void Tree::print(std::string header) {

    std::cout << header << std::endl;
    print();
}

void Tree::setBranch(Node* e1, Node* e2, double x) {
    std::pair<Node*,Node*> key;
    initializeBranchLengthKey(key, e1, e2);
    BranchLengths::iterator it = branchLengths.find(key);
    if (it != branchLengths.end())
        {
        it->second = x;
        }
    else
        {
        branchLengths.insert( std::make_pair(key, x) );
        }
}

void Tree::showNode(Node* p, int indent) {

    if (p != nullptr)
        {
        std::set<Node*>& pNeighbors = p->getNeighbors();
        for (int i=0; i<indent; i++)
            std::cout << " ";
        std::cout << p->getIndex() << " [" << p << "] ( ";
        for (Node* n : pNeighbors)
            {
            if (n == p->getAncestor())
                std::cout << "a_";
            std::cout << n->getIndex() << " ";
            }
        std::cout << ") ";
        std::cout << p->getName() << " " << p->getIsTip() << " ";
        std::cout << std::fixed << std::setprecision(6);
        if (p != root)
            std::cout << getBranchLength(p, p->getAncestor()) << " ";
        else
            std::cout << "--- ";
        if (p == root)
            std::cout << "<- Root ";
        std::cout << std::endl;

        for (Node* n : pNeighbors)
            {
            if (n != p->getAncestor())
                showNode(n, indent + 3);
            }
            
        }
}

void Tree::writeTree(Node* p, std::stringstream& strm) {

    if (p != nullptr)
        {
        std::set<Node*>& pNeighbors = p->getNeighbors();
        if (p->getIsTip() == false)
            strm << "(";
        else
            {
            std::string tName = p->getName();
            std::replace( tName.begin(), tName.end(), ' ', '_');
            strm << tName;
            }
        bool firstNode = true;
        for (Node* n : pNeighbors)
            {
            if (n != p->getAncestor())
                {
                if (firstNode == false)
                    strm << ",";
                writeTree(n, strm);
                strm << ":" << getBranchLength(p, n);
                firstNode = false;
                }
            }
        if (p->getIsTip() == false)
            strm << ")";

        }
}
