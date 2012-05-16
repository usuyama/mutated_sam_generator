@region = "chr21:17438150-17438299"
@chr = "chr21"
LEFT_POS = 17438150
REF  = "TTAATTCCATCAGGATCCTTAATTATCCTTTGCTACATAATGTAATATATGCAGTTTCTAGGGATCAGGATGTGGAAAATTTTGGAGGACCATCTATCTGTCACCTAAGAGGAAATAACTCAGAAGAGGCTGTTGAAACCAAAAAGCAA"
DNA = ["A", "T", "G", "C"]
ERROR_RATE = 0.01
READ_SIZE=75
a = ARGV[0].to_f
N = ARGV[1].to_i
@fn = [0.5, 0.5]
@ft = [0.5, 0.5-a, a/2.0, a/2.0]
p @fn
p @ft
class Variant
	attr_reader :type, :ref_pos, :length, :str, :seq
	attr_accessor :hap_pos
	def initialize(pos, str)
		@str = str
		@ref_pos = pos
		@hap_pos = pos
		if /(\w)=>(\w)/ =~ str
			@type = :snp
			@seq = $2
			@length = 1
		elsif /\+(\w+)/ =~ str
			@type = :ins
			@seq = $1
			@length = $1.length
		elsif /-(\w+)/ =~ str
			@type = :del
			@seq = $1
			@length = $1.length
		end
	end
	def inspect
		"#{@ref_pos}, #{@hap_pos} " + @str
	end
end

class Haplotype
	attr_reader :variants, :seq, :punc
	def initialize(vars=[])
		@variants = vars.sort {|a,b| a.ref_pos <=> b.ref_pos }
		gen_seq
		gen_punc
	end

	def inspect
		@seq
	end

	def gen_punc
		@punc = []
		for v in @variants
			if v.type == :snp
				next	
			elsif v.type == :ins
				@punc.push([:is, v.hap_pos])
				@punc.push([:ie, v.hap_pos+v.length])
			elsif v.type == :del
				@punc.push([:ds, v.hap_pos])
				@punc.push([:de, v.hap_pos+v.length])
			end
		end
	end

	def change_hap_pos(start, length)
		start.upto(@variants.size-1) do |i|
			@variants[i].hap_pos += length
		end
	end

	def gen_seq
		@seq = REF.dup
		i = 0
		for v in @variants
			i += 1
			if v.type == :snp
				@seq[v.hap_pos-1] = v.seq
			elsif v.type == :ins
				@seq = @seq[0..v.hap_pos-1] + v.seq + @seq[v.hap_pos..@seq.size-1]
				change_hap_pos(i, v.length)
			elsif v.type == :del
				@seq = @seq[0..v.hap_pos-1] + @seq[v.hap_pos+v.length..@seq.size-1]
				change_hap_pos(i, -v.length)
			end
		end	
	end

	def seek_next(pos, rest)
		left = nil
		right = nil
		old = nil
		for x in @punc
			if x.last > pos
				right = x
				if !old.nil?
					left = old
				end
				break
			end					
			old = x
		end
		if left.nil?
			if !right.nil?
				l = right.last - pos
				l = (l < rest) ? l : rest
				return [l, "#{l}M"]	
			else
				[rest, "#{rest}M"]	
			end
		end		
		if right.nil?
			[rest, "#{rest}M"]
		else
			l = right.last - pos
			l = l > rest ? rest : l
			if  left.first == :ie || left.first == :de
				[l, "#{l}M"]			
			elsif left.first == :is
				[l, "#{l}I"]
			else
				l = right.last - pos
				[l, "#{l}D"]
			end
		end
	end

	def gen_cigar(start)
		cp = start
		rest = READ_SIZE
		cigar = ""
		while(rest > 0)
			x = seek_next(cp, rest)
			cigar += x.last
			if x.last[-1,1]!="D"
				rest -= x.first
			end
			cp += x.first
		end
		cigar
	end


	def gen_read
		pos = rand(REF.length-READ_SIZE-1)
		read = @seq[pos, READ_SIZE]
		cigar = gen_cigar(pos)
		0.upto(read.size-1) do |i|
			if(rand() < ERROR_RATE)
				alt = read[i,1]
				while(alt==read[i,1])
					alt = DNA[rand(4)]
				end
				read[i] = alt
			end
		end
		[pos+1, read, cigar, "C"*READ_SIZE]
	end
end
@h1 = Haplotype.new
p @h1
@h2 = Haplotype.new([Variant.new(83,"-GGA")])
@h3 = Haplotype.new([Variant.new(86,"-GGA")])
@h4 = Haplotype.new([Variant.new(83,"T=>A")])
@hn = [	@h1,@h4]
@ht = [	@h1,@h4, @h2,@h3]



p @hn
p @ht
@normal_reads = []
N.times do 
	r = rand()
	i = 0
	th = 0.0
	1.upto(@fn.size()) do
		th += @fn[i]
		if th > r
			break
		else	
			i += 1
		end
	end
	@normal_reads.push(@hn[i].gen_read())
end

@tumor_reads = []
N.times do 
	r = rand()
	i = 0
	th = 0.0
	1.upto(@ft.size()) do
		th += @ft[i]
		if th > r
			break
		else	
			i += 1
		end
	end
	@tumor_reads.push(@ht[i].gen_read())
end
@header
File.open(File.expand_path(File.dirname($0)) + "/header", "r") do |f|
	@header = f.read()
end

File.open("out_normal", "w") do |f|
	f.write(@header)
	for r in @normal_reads.sort {|a,b| a<=>b }
		f.write("USU_SIM\t0\t#{@chr}\t#{LEFT_POS+r[0]}\t100\t#{r[2]}\t*\t0\t0\t#{r[1]}\t#{r[3]}\n")
	end
end

File.open("out_tumor", "w") do |f|
	f.write(@header)
	for r in @tumor_reads.sort {|a,b| a<=> b}
		f.write("USU_SIM\t0\t#{@chr}\t#{LEFT_POS+r[0]}\t100\t#{r[2]}\t*\t0\t0\t#{r[1]}\t#{r[3]}\n")
	end
end

