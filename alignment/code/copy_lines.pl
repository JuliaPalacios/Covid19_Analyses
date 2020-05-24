#!/usr/bin/perl

# From https://unix.stackexchange.com/a/506226
# usage: thisscript linenumberslist.txt contentsfile

unless (open(IN, $ARGV[0])) {
        die "Can't open list of line numbers file '$ARGV[0]'\n";
}
my %linenumbers = ();
while (<IN>) {
        chomp;
        $linenumbers{$_} = 1;
}

unless (open(IN, $ARGV[1])) {
        die "Can't open contents file '$ARGV[1]'\n";
}
$. = 0;
while (<IN>) {
        print if defined $linenumbers{$.};
}

exit;