#ifndef Node_hpp
#define Node_hpp

#include <set>
#include <string>
class RandomVariable;


class Node {

    public:
                            Node(void);
        void                addNeighbor(Node* p);
        Node*               getAncestor(void) { return ancestor; }
        std::vector<Node*>& getDescendants(void);
        int                 getIndex(void) { return index; }
        bool                getIsTip(void) { return isTip; }
        std::string         getName(void) { return name; }
        std::set<Node*>&    getNeighbors(void) { return neighbors; }
        int                 getOffset(void) { return offset; }
        void                removeNeighbor(Node* p);
        void                removeAllNeighbors(void) { neighbors.clear(); }
        void                setAncestor(Node* p) { ancestor = p; }
        void                setIndex(int x) { index = x; }
        void                setIsTip(bool tf) { isTip = tf; }
        void                setName(std::string s) { name = s; }
        void                setOffset(int x) { offset = x; }
    
    private:
        Node*               ancestor;
        double              branchLength;
        std::vector<Node*>  descendantsVector;
        int                 index;
        bool                isTip;
        std::string         name;
        std::set<Node*>     neighbors;
        int                 offset;
};

#endif
