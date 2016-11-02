#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

char ch[2000000];
string tmp;

int main()
{
    freopen("1/testnet3/debug.log", "r", stdin);
//    freopen("BLOCKS.LOG", "w", stdout);
    FILE *blk = fopen("Blocks.log", "w");
    FILE *frt = fopen("Fruits.log", "w");
    while (gets(ch)) {
        tmp = string(ch);
        if (tmp.find("Block found") != tmp.npos) fprintf(blk, "%s\n", tmp.c_str());
        else if (tmp.find("Fruit found") != tmp.npos) fprintf(frt, "%s\n", tmp.c_str());
    }
}
