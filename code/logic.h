#include <string>
#include <vector>

using namespace std;

class Prop {
public:
    enum Connective {
        AND,
        OR,
        IMP,
        TOP,
        BOT,
        ATOM
    };
    Connective conn;
    std::string name;
    Prop *left;
    Prop *right;
    Prop(Connective c, Prop* l = nullptr, Prop* r = nullptr)
        : conn(c), left(l), right(r) {}
    Prop(std::string n)  // atomic proposition
        : conn(ATOM), name(std::move(n)), left(nullptr), right(nullptr) {}
};

string toString(const Prop* p);
bool rightInversion(vector<Prop*> &db, vector<Prop*> &stoup, const Prop *goal);
bool leftInversion(vector<Prop*> &db, vector<Prop*> &stoup, const Prop *goal);
bool searchRight(vector<Prop*> &db, const Prop *goal);
bool searchLeft(vector<Prop*> &db, const Prop *goal);
bool prove(const Prop *goal);