package Bio::Align::Subset;

use 5.006;
use strict;
use warnings;

=head1 NAME

Bio::Align::Subset - A BioPerl module to generate new alignments as subset from larger alignments

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Given an array of codon positions and an alignment, the function alignSubset
returns a new alignment with the codons at those positions from the original
alignment.

Usage example:
    
    use Bio::Align::Subset;
    
    # Align file name
    my $inAlignFile = "align.fasta";
    # Load the file content in a file handle (FH)
    open(FH, $inAlignFile);
    # Align data
    my @alignData = <FH>;
    
    # Build the codon sets
    my @set1 = ();
    my @set2 = ();
    my @set3 = ();
    
    # Store the set references
    my $setRef1 = \@set1;
    my $setRef2 = \@set2;
    my $setRef3 = \@set3;
    
    # Construct a new object
    my $subAlign = Bio::Align::Subset->new();
    # Pass set references by arguments
    $subAlign->subsetsReferences = [$setRef1, $setRef2, $setRef3];
    # Pass the original alignment data
    $subAlign->originalAlignFH = @alignData;
    # Generates the new alignments as strings
    my @subSets = $subAlign->launch();
    
    # Store each subset in a file
    for(my $i=0;$i<=$#subSets;$i++){
        open("FOUT".$i,">file".$i.".fasta");
        print FOUT, $subSets[$i];
        fclose(FOUT);
    }
    
    

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 alignSubset

=cut

sub alignSubset {
    
    
    
}



=head1 AUTHOR

Hector Valverde and Juan Carlos Aledo, C<< <hvalverdeuma.es> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-align-subset at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Align-Subset>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Align::Subset


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Align-Subset>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Align-Subset>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Align-Subset>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Align-Subset/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2012 Hector Valverde and Juan Carlos Aledo.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of Bio::Align::Subset