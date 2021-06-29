
function make_fqparser( filename; forcegzip=false )
   fopen = open( filename, "r" )
   if isgzipped( filename ) || forcegzip
      to_open = ZlibInflateInputStream( fopen, reset_on_end=true )
   else
      to_open = BufferedInputStream( fopen )
   end
   FASTQ.Reader( to_open, fill_ambiguous=DNA_A )
end

# modified for Bio v0.2 with tryread_bool!
@inline function read_chunk!( chunk, parser )
   i = 1
   while i <= length(chunk) && !eof(parser)
      read!( parser, chunk[i] )
      i += 1
   end
   while i <= length(chunk)
      pop!(chunk) # clean up if we are at the end
   end
   parser
end

function allocate_chunk( parser; size=10000 )
  chunk = Vector{eltype(parser)}( size )
  for i in 1:length(chunk)
     chunk[i] = eltype(parser)()
  end
  chunk
end

function allocate_fastq_records( size::Int=10000 )
   chunk = Vector{FASTQRecord}( undef, size )
   for i in 1:length(chunk)
      chunk[i] = FASTQRecord()
   end
   chunk
end

function process_reads!( parser, param::AlignParam, lib::GraphLib, quant::GraphLibQuant,
                         multi::MultiMapping{SGAlignSingle}, mod::B;
                         bufsize=150, sam=false, qualoffset=33 ) where B <: BiasModel

#   println("||||$(param)||||||||$(quant)||||")

   reads  = allocate_fastq_records( bufsize )
   mean_readlen = 0.0
   total        = 0
   mapped       = 0
   if sam
      stdbuf = BufferedOutputStream( stdout )
      write_sam_header( stdbuf, lib )
   end

   fqfile = open("trash.fastq", "w")
   fqwriter = FASTQ.Writer(fqfile)

   while length(reads) > 0
      read_chunk!( reads, parser )
      total += length(reads)
      println(stderr, "$(length(reads)) reads/$total............")
      @inbounds for i in 1:length(reads)
#= Lots of these 'reads' objects have an 'EMPTY SEQUENCE' and the ones that do not
   appear to have a sequence that does not match the .sequence data.  Here's one.

Whippet.FASTQRecord(< EMPTY SEQUENCE >, UInt8[], FASTX.FASTQ.Record:
   identifier: HISEQ:587:C6RMKANXX:1:1307:3191:12232/1
  description: <missing>
     sequence: CTCAGAGTGAAAGGCTGGAAAACAAATTTCCAAGCAAATGGTCTGAAGAAACAAGCTGGAGTAGCCATTCTAATATCGAATAAAATCGACTTCCAACCCAAAGTTATCAAAAAAGACAAGGAGGGA
      quality: UInt8[...])
=#
#         println(stderr, ">>>$(reads[i])")
#         println(stderr, "***$(reads[i].raw)")
#	 println(stderr, "^^^$(first(reads[i].raw.sequence))^^^$(last(reads[i].raw.sequence))")
         fill!( reads[i], qualoffset )
#= The 'fill' this refers to appears to be that it fills in the EMPTY SEQUENCE with
   the sequence in the '.sequence' field.
=#
#	 println(stderr, "<<<$(reads[i])")
         align = ungapped_align( param, lib, reads[i] )
         if !isnull( align )
            descString = ""
	    for k in 1:length(align.value)
               for m in 1:length(align.value[k].path)
                  descString *= lib.names[align.value[k].path[m].gene] * ":";
                  descString *= string(align.value[k].path[m].node) * "(";
                  descString *= string(align.value[k].path[m].score.matches) * ",";
                  descString *= string(align.value[k].path[m].score.mismatches) * ",";
                  descString *= string(align.value[k].path[m].score.mistolerance) * ")/";
               end
            end
            ## Writes a FASTQ file for any reads that match the given gene.
            if lib.names[first(first(align.value).path).gene] == "ENSMUSG00000040653"
               write(fqwriter,
                     FASTQ.Record(identifier(reads[i].raw),
                                  descString,
                                  sequence(reads[i].raw),
                                  quality(reads[i].raw); offset=33));
            end
            biasval = count!( mod, reads[i].sequence )
            if length( align.value ) > 1
               push!( multi, align.value, biasval, quant, lib )
               sam && write_sam( stdbuf, reads[i], align.value, lib, qualoffset=qualoffset )
            else
               count!( quant, align.value[1], biasval )
               sam && write_sam( stdbuf, reads[i], align.value[1], lib, qualoffset=qualoffset )
            end
            mapped += 1
            @fastmath mean_readlen += (length(reads[i].sequence) - mean_readlen) / mapped
         end
      end
      if total % 100000 == 0
         GC.gc()
      end
   close(fqfile)
   end # end while
   if sam
      close(stdbuf)
   end
   mapped,total,mean_readlen
end


function process_paired_reads!( fwd_parser, rev_parser, param::AlignParam,
                                lib::GraphLib, quant::GraphLibQuant,
                                multi::MultiMapping{SGAlignPaired}, mod::B;
                                bufsize=50, sam=false, qualoffset=33 ) where B <: BiasModel

   fwd_reads  = allocate_fastq_records( bufsize )
   rev_reads  = allocate_fastq_records( bufsize )
   mean_readlen = 0.0
   total        = 0
   mapped       = 0
   if sam
      stdbuf = BufferedOutputStream( stdout )
      write_sam_header( stdbuf, lib )
   end
   while length(fwd_reads) > 0 && length(rev_reads) > 0
      read_chunk!( fwd_reads, fwd_parser )
      read_chunk!( rev_reads, rev_parser )
      total += length(fwd_reads)
      @inbounds for i in 1:length(fwd_reads)
         fill!( fwd_reads[i], qualoffset )
         fill!( rev_reads[i], qualoffset )
         fwd_aln,rev_aln = ungapped_align( param, lib, fwd_reads[i], rev_reads[i] )
         if !isnull( fwd_aln ) && !isnull( rev_aln )
            biasval = count!( mod, fwd_reads[i].sequence, rev_reads[i].sequence )
            if length( fwd_aln.value ) > 1
               push!( multi, fwd_aln.value, rev_aln.value, biasval, quant, lib )
               sam && write_sam( stdbuf, fwd_reads[i], rev_reads[i], fwd_aln.value, rev_aln.value, lib,
                                 paired=true, is_pair_rc=param.is_pair_rc, qualoffset=qualoffset )
            else
               count!( quant, fwd_aln.value[1], rev_aln.value[1], biasval )
               sam && write_sam( stdbuf, fwd_reads[i], fwd_aln.value[1], lib,
                                 paired=true, fwd_mate=true, is_pair_rc=param.is_pair_rc,
                                 qualoffset=qualoffset )
               sam && write_sam( stdbuf, rev_reads[i], rev_aln.value[1], lib,
                                 paired=true, fwd_mate=false, is_pair_rc=param.is_pair_rc,
                                 qualoffset=qualoffset )
            end
            mapped += 1
            @fastmath mean_readlen += (length(fwd_reads[i].sequence) - mean_readlen) / mapped
         end
      end
   end # end while
   if sam
      close(stdbuf)
   end
   mapped,total,mean_readlen
end
