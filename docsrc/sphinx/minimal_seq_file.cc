#include "reprimand/unitconv.h"
#include "reprimand/star_sequence.h"
#include "reprimand/star_seq_file.h"


using namespace EOS_Toolkit;

int main() {
    auto u = units::geom_solar();

    auto seq1 = load_star_branch("example1.tovseq.h5", u);    
    save_star_branch("copy1.tovseq.h5", seq1);

    auto seq2 = load_star_seq("example2.tovseq.h5", u);    
    save_star_seq("copy2.tovseq.h5", seq2);

}
