#ifndef Tree_hpp
#define Tree_hpp

#include <map>
#include <sstream>
#include <string>
#include <vector>
class Node;
typedef std::map<std::pair<Node*,Node*>,double> BranchLengths;


class Tree {

    public:
                                            Tree(void) = delete;
                                            Tree(const Tree& t); //copy constructor
                                            Tree(std::string newick); //tree from newick string
                                           ~Tree(void);
        Tree&                               operator=(const Tree& t);
        void                                checkBranchLengthsNeg(void);
        double                              getBranchLength(Node* e1, Node* e2);
        std::vector<Node*>&                 getDownPassSequence(void) { return downPassSequence; }
        std::string                         getNewickString(void);
        int                                 getNumNodes(void) { return (int)nodes.size(); }
        int                                 getNumTaxa(void) { return numTaxa; }
        Node*                               getRoot(void) { return root; }
        void                                initializeDownPassSequence(void);
        void                                print(void);
        void                                print(std::string header);
        void                                setBranch(Node* e1, Node* e2, double x);
        
    private:
        Node*                               addNode(void);
        void                                clone(const Tree& t);
        void                                deleteNodes(void);
        void                                initializeBranchLengthKey(std::pair<Node*,Node*>& key, Node* e1, Node* e2);
        std::vector<std::string>            parseNewickString(std::string);
        void                                passDown(Node* p, Node* from);
        void                                showNode(Node* p, int indent);
        void                                writeTree(Node* p, std::stringstream& strm);
        BranchLengths                       branchLengths;
        std::vector<Node*>                  downPassSequence;
        Node*                               freeNode;
        std::vector<Node*>                  nodes;
        int                                 numTaxa;
        Node*                               root;
        std::vector<Node*>                  tips;
};

#endif
