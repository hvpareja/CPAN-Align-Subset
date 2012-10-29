package Bio::Align::Subset;

use 5.006;
use strict;
use warnings;
use Carp;

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

=cut

# Body ########################################################################
###############################################################################
###############################################################################
###############################################################################

=head1 CONSTRUCTOR

=head2 Bio::Align::Subset->new

=cut
#
# The constructor of the class
#
sub new {
    
    my ($class, %arg) = @_;
    my $self = bless {
        _inputfile => $arg{inputfile} || croak("No file specified"),
        _subsets   => $arg{subsets} || croak("No subsets defined")
    }, $class;
    $class->_incr_count();
    return $self;
    
}

###############################################################################
# Tracking methods
###############################################################################
{
    my $_count = 0;
    sub get_count{
        $_count;
    }
    sub _incr_count{
        ++$_count;
    }
    sub _decr_count{
        --$_count;
    }
}
###############################################################################

###############################################################################
# Auxiliary methods
###############################################################################
{
    
    my $_sequence_length = 0;
    
    #
    # Set the sequence length of the whole alignment
    #
    sub _set_sequence_length{
        $_sequence_length = $_[0];
    }
    
    #
    # Check if a the length of a given sequence match with the length of
    # the whole alignment.
    #
    sub _check_sequence_length{
        my $tested_sequence_length = $_[0];
        $tested_sequence_length == $_sequence_length ? return 1 : return 0;
    }
    
    #
    # Verifies the integrity of a given sequence
    #
    sub _verify_chain{
        
        my $sequence = $_[0];
        my $seq_length = length($sequence);
        
        # 1. The chain must be a DNA sequence
        _isdna($sequence) ? 1 : croak("The sequence is not a dna sequence");
        
        # 2. Also, the input file must be wrapped (non untermitated codons)
        $seq_length % 3 != 0 ? 1 : croak("The sequence length is not multiple of 3");
        
        # 3. Finally, all the sequences must be equal. But if $_sequence_length
        # has not been updated, it takes the value of the length of this sequence.
        if($_sequence_length == 0){
            _set_sequence_length($seq_length);
        }else{
            _check_sequence_length($seq_length) ? 1 : croak("A sequence length does not match with le length of the whole alignment");
        }
        
        return 1;
        
    }
    
    #
    # Verifies if a given string is a DNA sequence
    #
    sub _isdna{
        
        my $sequence = uc($_[0]);
        $sequence =~/^[ACGT]+$/i ? return 1: return 0;
    }
    
    
    
}
###############################################################################


###############################################################################
# Accessor Methods
###############################################################################
# This kind of method is called Accesor
# Method. It returns the value of a key
# and avoid the direct acces to the inner
# value of $obj->{_inputfile}.
###############################################################################
sub get_inputfile { $_[0] -> {_inputfile} } 
sub get_subsets { $_[0] -> {_subsets} }
###############################################################################


###############################################################################
# Mutator Methods
###############################################################################
sub set_inputfile { my ($self, $inputfile) = @_;
                    $self-> {_inputfile} = $inputfile if $inputfile;
                  }
sub set_subsets   { my ($self, $subsets) = @_;
                    $self-> {_subsets} = $subsets if $subsets;
                  }
###############################################################################



# Footer ######################################################################
###############################################################################
###############################################################################
###############################################################################

=head1 AUTHOR

Hector Valverde and Juan Carlos Aledo, C<< <hvalverde@uma.es> >>

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
