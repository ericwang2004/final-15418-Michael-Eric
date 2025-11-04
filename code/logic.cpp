#include "logic.h"
#include <iostream>

using namespace std;

string toString(const Prop *p) {
    switch (p->conn) {
        case Prop::ATOM:
            return p->name;
        case Prop::AND:
            return "(" + toString(p->left) + " /\\ " + toString(p->right) + ")";
        case Prop::OR:
            return "(" + toString(p->left) + " \\/ " + toString(p->right) + ")";
        case Prop::IMP:
            return "(" + toString(p->left) + " => " + toString(p->right) + ")";
        case Prop::TOP:
            return "T";
        case Prop::BOT:
            return "F";
        default:
            exit(EXIT_FAILURE);
    }
}

bool searchAtomic(vector<Prop*> &db, string name) {
    for (size_t i = 0; i < db.size(); i++) {
        Prop *P = db[i];
        if (P->conn == Prop::ATOM && P->name == name) {
            return true;
        }
    }
    return false;
}

bool rightInversion(vector<Prop*> &db, vector<Prop*> &stoup, const Prop *goal) {
    if (!goal) {
        exit(EXIT_FAILURE);
    }
    switch (goal->conn) {
        case Prop::AND:
            return rightInversion(db, stoup, goal->left) && rightInversion(db, stoup, goal->right);
        case Prop::TOP:
            return true;
        case Prop::IMP: {
            stoup.push_back(goal->left);
            return rightInversion(db, stoup, goal->right);
        }
        default:
            return leftInversion(db, stoup, goal);
    }
}

bool leftInversion(vector<Prop*> &db, vector<Prop*> &stoup, const Prop *goal) {
    if (stoup.empty()) {
        return searchRight(db, goal);
    }
    Prop *focus = stoup.back();
    stoup.pop_back();
    switch (focus->conn) {
        case Prop::AND: {
            stoup.push_back(focus->left);
            stoup.push_back(focus->right);
            return leftInversion(db, stoup, goal);
        }
        case Prop::TOP: 
            return leftInversion(db, stoup, goal);
        case Prop::OR: {
            vector<Prop*> stoup1 = stoup;
            stoup1.push_back(focus->left);
            vector<Prop*> stoup2 = stoup;
            stoup2.push_back(focus->right);
            return leftInversion(db, stoup1, goal) && leftInversion(db, stoup2, goal);
        }
        case Prop::BOT:
            return true;
        case Prop::IMP: {
            Prop *ante = focus->left;
            switch (ante->conn) {
                case Prop::TOP: {
                    stoup.push_back(focus->right);
                    return leftInversion(db, stoup, goal);
                }
                case Prop::BOT:
                    return leftInversion(db, stoup, goal);
                case Prop::AND: {
                    Prop *a = ante->left;
                    Prop *b = ante->right;
                    Prop *c = focus->right;
                    Prop *curry = new Prop(Prop::IMP, a, new Prop(Prop::IMP, b, c));
                    stoup.push_back(curry);
                    return leftInversion(db, stoup, goal);
                }
                case Prop::OR: {
                    Prop *a = ante->left;
                    Prop *b = ante->right;
                    Prop *c = focus->right;
                    Prop *ac = new Prop(Prop::IMP, a, c);
                    Prop *bc = new Prop(Prop::IMP, b, c);
                    stoup.push_back(bc);
                    stoup.push_back(ac);
                    return leftInversion(db, stoup, goal);
                }
                default: {
                    db.push_back(focus);
                    return leftInversion(db, stoup, goal);
                }
            }
        }
        default: {
            db.push_back(focus);
            return leftInversion(db, stoup, goal);
        }
    }
}

bool searchRight(vector<Prop*> &db, const Prop *goal) {
    switch (goal->conn) {
        case Prop::OR: {
            vector<Prop*> empty;
            if (rightInversion(db, empty, goal->left) || rightInversion(db, empty, goal->right)) {
                return true;
            }
            return searchLeft(db, goal);
        }
        case Prop::ATOM: {
            if (searchAtomic(db, goal->name)) {
                return true;
            }
            return searchLeft(db, goal);
        }
        default:
            return searchLeft(db, goal);
    }
}

bool searchLeft(vector<Prop*> &db, const Prop *goal) {
    if (db.empty()) {
        return false;
    }
    vector<Prop*> db_copy = db;
    for (size_t i = 0; i < db.size(); i++) {
        Prop *focus = db[i];
        db.erase(db.begin() + i);
        
    }
}

bool prove(const Prop *goal) {
    vector<Prop*> db;
    vector<Prop*> stoup;
    return rightInversion(db, stoup, goal);
}

int main() {
    Prop P("P");
    // And(T, Or(T, R))
    // Prop *test1 = new Prop(Prop::AND, new Prop(Prop::TOP), new Prop(Prop::AND, new Prop(Prop::TOP), new Prop(Prop::TOP)));
    Prop *test2 = new Prop(Prop::IMP, new Prop(Prop::BOT), &P);
    std::cout << toString(test2) << ": " << prove(test2) << "\n";
    return 0;
}