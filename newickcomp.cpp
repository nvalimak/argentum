/**
 * Helper script to parse Newick trees
 *
 * Usage:
 *         ./main [options] | ./newickcomp <SCRM>
 *
 * Input: Predicted Newick trees from 'main' (as standard input).
 *        <SCRM> SCRM output file name.
 *
 * Example:
 *        ./main --scrm --input trees100.txt --verbose --newick --nrow 1000 | ./newickcomp trees100.txt
 */
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <cassert>
#include <cstdlib>
using namespace std;

int atoi_min(char const *value, int min)
{
    std::istringstream iss(value);
    int i;
    char c;
    if (!(iss >> i) || iss.get(c))
    {
        cerr << "error: value " << value << " must be of type <int>, and greater than or equal to " << min << endl;
        std::exit(1);
    }

    if (i < min)
    {
        cerr << "error: value " << value << " must be greater than or equal to " << min << endl;
        std::exit(1);
    }
    return i;
}

class NewickTree
{
public:
    class Node
    {
    public:
        Node()
            : ch(), parent(0), label(-1), leaf(false), size(0)
        { }
        Node(Node *p)
            : ch(), parent(p), label(-1), leaf(false), size(0)
        {
            p->ch.insert(this); // Insert as child to *p
        }
        
        ~Node()
        { for (set<Node *>::iterator it = ch.begin(); it != ch.end(); ++it) delete *it; }
        set<Node *> ch; // Children
        Node *parent; 
        int label; // Leaf label
        bool leaf; // Is leaf?
        unsigned size; // Total number of subtree leaves
    };

    /**
     * Update the distance values of 'src_leaf' by traversing whole tree
     */
    void distanceByTraversal(Node *pn, Node *src, Node *src_leaf, unsigned maxval, vector<unsigned> &cd)
    {
        maxval = max(maxval, pn->size);
        if (pn->leaf)
        {
            unsigned i = pn->label;
            unsigned j = src_leaf->label;
            assert (i < root->size && j < root->size);
            cd[i + root->size * j] = maxval;
            return;
        }

        for (set<Node *>::iterator it = pn->ch.begin(); it != pn->ch.end(); ++it)
            if (*it != src) // Assert: src has already been processed
                distanceByTraversal(*it, pn, src_leaf, maxval, cd);
        if (pn->parent != 0 && pn->parent != src)
            distanceByTraversal(pn->parent, pn, src_leaf, maxval, cd);
    }
    /**
     * Update distance values over all leaf node pairs
     */
    void subtreeDistance(vector<unsigned> &cd)
    {
        assert (cd.size() == root->size*root->size);
        for (unsigned i = 0; i < cd.size(); ++i)
            cd[i] = 0;
        for (unsigned i = 0; i < leaves.size(); ++i)
            distanceByTraversal(leaves[i]->parent, leaves[i], leaves[i], 1, cd);
    }

    /* FIXME Not needed? 
    Node * commonAncestor(Node *pn, Node *rn)
    {
        set<Node *> pn_parents;
        while (pn->parent != 0)
        {
            pn_parents.insert(pn->parent);
            pn = pn->parent;
        }
        do 
            rn = rn->parent;
        while (pn_parents.count(rn) == 0);
        return rn;
        }*/

    // Create a new tree from the given string
    NewickTree(string const &i)
        : root(0), valid(-1)
    {
        parse(i);
        updateSizes(root);
        assert (root->size == leaves.size());
    }

    // Bottom-up cascade of subtree size (n:o leaves)
    int updateSizes(Node *pn)
    {
        if (pn->leaf)
        {
            pn->size = 1;
            return 1;
        }
        pn->size = 0;
        for (set<Node *>::iterator it = pn->ch.begin(); it != pn->ch.end(); ++it)
            pn->size += updateSizes(*it);
        return pn->size;
    }

    // Simple tree parser
    size_t parse(string const &s, size_t i, Node *p)
    {
        while (i != string::npos && i < s.size()-1)
        {
            if (s[i]=='(')
            {
                Node *j = 0;
                if (p == 0)
                {
                    j = new Node();
                    root = j;
                }
                else
                    j = new Node(p);
                i = parse(s, i+1, j);
                continue;
            }
            else if (s[i] == ',')
            {
                i+=1;
                continue;
            }
            else if (s[i] == ')')
            {
                i = s.find_first_of(",)(;",i+1);
                return i;
            }
            else
            {
                Node *j = new Node(p);
                j->label = ::atoi_min(s.substr(i, s.find_first_of(":",i)-i).c_str(), 1)-1;
                j->leaf = true;
                leaves.push_back(j);
                i = s.find_first_of("),",i);
            }
        }
        return i;
    }
    void parse(string const &s)
    {
        if (s[0] != '[')
            cerr << "error on row: " << s << endl;
        assert(s[0] == '[');
        size_t i = s.find_first_of("]");
        valid = ::atoi_min(s.substr(1,i-1).c_str(), 0);
        i +=1;
        i = parse(s, i, 0);
        assert(s[i] == ';');
    }

    /** 
     * Debuging output (Graphwiz DOT format)
     */
    unsigned outputDOT(Node *pn, unsigned id)
    {
        unsigned oid = id;
        if (pn->leaf)
            return id;
        for (set<Node *>::iterator it = pn->ch.begin(); it != pn->ch.end(); ++it)
        {
            id ++;
            cout << " n" << oid << " -> n" << id << endl;
            cout << " n" << id << " [label=\"";
            if ((*it)->leaf)
                cout << (*it)->label;
            else
                cout << id;
            cout << "\"]" << endl;
            id = outputDOT(*it, id);
        }
        return id;
    }
    /** 
     * Debuging output (Graphwiz DOT format)
     */
    void outputDOT()
    {
        cout << "digraph G {" << endl;
        outputDOT(root, 100);
        cout << "}" << endl;
    }

    Node *root;
    vector<Node *> leaves;
    int valid; // Tree is valid for this many sites
};

    
int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        cerr << "usage: ./main [options] | " << argv[0] << " <SCRM-input-file>" << endl;
        return 1;
    }

    ifstream scrmf(argv[1]);
    string row;
    if (!std::getline(scrmf, row).good())
    {
        cerr << "error: unable to read the input file " << argv[1] << endl;
        return 1;
    }
    while (row[0] != '[') // Find first tree
        std::getline(scrmf, row);
    NewickTree * tree = new NewickTree(row);
    unsigned nleaves = tree->root->size;
    vector<unsigned> cd_orig(nleaves*nleaves, 0);
    tree->subtreeDistance(cd_orig);
    unsigned site = 0;
    unsigned n = 0;
    while (scrmf.good())
    {
        if (site >= tree->valid)
        {
            std::getline(scrmf, row);
            if (!scrmf.good() || row[0] != '[')
            {
                cerr << "EOF from scrm input after n = " << n << " sites" << endl;
                return 0;
            }
            site = 0;
            delete tree;
            // update tree & distances
            tree = new NewickTree(row);
            tree->subtreeDistance(cd_orig);
            assert (tree->root->size == nleaves);
        }

        if (!std::getline(cin, row).good())
            break;
        if (row.empty())
            break;
        NewickTree *predt = new NewickTree(row);
        assert (predt->root->size == tree->root->size);
        vector<unsigned> cd(nleaves*nleaves, 0);
        predt->subtreeDistance(cd);

        unsigned i = 0;
        unsigned j = 1;
        // Output <site> <leaf id> <leaf id> <dist. in SCRM> <dist. in stdin>
        cout << n << '\t' << i+1 << '\t' << j+1 << '\t' << cd_orig[i + nleaves * j] << '\t' << cd[i + nleaves * j] << endl;
        
        delete predt; predt = 0;
        site ++;
        n ++;
    }
    delete tree;
    return 0;
}
