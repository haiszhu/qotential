#!/usr/bin/env perl
use strict;
use warnings;

# LFortran 0.64 can print a complete LLVM module and then append its verifier
# diagnostic.  Keep only the module, and remove the redundant integer type on
# literal operands such as `add i64 %idx, i32 1` -> `add i64 %idx, 1`.
# LFortran 0.65 also emits the non-semantic `nocreateundeforpoison` and
# `captures(none)` optimizer attributes, which Apple clang 17 cannot parse;
# omit them for host compilation.
# The rewrite assumes LLVM's textual form of one instruction per line.  If a
# future LFortran emits something else, fail loudly rather than hand clang a
# quietly corrupted module.
my $nline = 0;
my $seen_define = 0;
while (my $line = <STDIN>) {
    last if $line =~ /^asr_to_llvm:/;
    $seen_define = 1 if $line =~ /^define /;
    if ($line =~ /\b(?:add|sub|mul|sdiv) i(?:32|64) /) {
        $line =~ s/, i(?:32|64) (-?[0-9]+)/, $1/g;
    }
    $line =~ s/\bnocreateundeforpoison\b\s*//g;
    $line =~ s/\s+captures\(none\)//g;
    print $line;
    $nline++;
}
die "fix-lfortran-llvm: no LLVM IR on stdin (did lfortran fail?)\n"
    if $nline == 0;
die "fix-lfortran-llvm: IR has no 'define' -- output looks truncated\n"
    unless $seen_define;
