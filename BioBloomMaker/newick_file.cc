
#include "newick_file.hh"
#include "newick.hh"
#include "tree.hh"
#include <sys/stat.h>
#include <unistd.h>

using namespace BiRC::treelib;

#include <sys/mman.h>
#include <fcntl.h>

// FIXME: There is no error handling here!
std::shared_ptr<Tree> BiRC::treelib::parse_newick_file(const char *fname) {
    const int filedes = open(fname, 0, 0); // FIXME: not exception safe
    struct stat file_stat;
    fstat(filedes, &file_stat);
    const size_t len = file_stat.st_size;
    const int prot = PROT_READ;
    const int flags = MAP_PRIVATE;
    const off_t off = 0;

    // FIXME: There is no error handling here!
    char *beg = static_cast<char *>(mmap(0, len, prot, flags, filedes, off));
    char *end = beg + len;

    std::shared_ptr<Tree> tree(parse_newick(beg, end));

    munmap(beg, len);
    close(filedes);

    return tree;
}

std::shared_ptr<Tree>
BiRC::treelib::parse_newick_file(const std::string &fname)
{
    return parse_newick_file(fname.c_str());
}
