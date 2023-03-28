ch_prova = Channel.fromList([["sample5_T1", "single_end"], ["/home/valentino/Program/gatkcaller-main/work/41/f1b5e17f1d759f3cb50fe588e98207/sample5_T1.md.bam", "/home/valentino/Program/gatkcaller-main/work/b3/6d855bec7a96dc44b5653f75cb4706/sample5_T1.dragstr_model.txt"]]
)
ch_prova
.map { id, the_list ->
        tuple( the_list )
    }
//.map {meta, bammodel -> tuple(meta, bammodel)}
//.map { tuple( it[0], it[1] ) }
//.map{
//    it ->
//    [it[0], it[1][0]]
//}
//.map { meta, meta2, bammodel -> [meta, bammodel] }
//.map { meta, bam -> (bam[1])}
//.filter{ it ->  (it.last())}
.view()

/*
    your_channel.map { id, the_list ->
        tuple( id, the_list.findAll { it.first() == "Type1" } )
    }
ch_alignment
  // For every element of this channel, convert it to a string, split in pieces separated by --, get the second part, then split by _3p and get the first part. Return a list with this as the first value, and then the original element as the second value. This part has to be customized depending on what part of the String you want to get as matching key
  .map { [it.toString().split("--")[1].split("_3p")[0],
          it] }.
  set { ch_alignment }
ch_clustered
  .map { [it.toString().split("--")[1].split("_3p")[0],
          it] }.
  set { ch_clustered }

ch_alignment
  // Combine according to a key that is the first value of every first element, which is a list according to what we did above
  .combine(ch_clustered, by: 0)
  // For every element of this channel, which consists of three values now, the matching key (id), the first element of the first channel, and the second, keep only the second and the third.
  .map { id, sam, fasta -> [sam, fasta] }
  // View the content of the channel, which consists of the last two values
  .view()

        ch_cram_variant_calling.branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }.set{ch_cram_variant_calling_status}

    .map {
        meta, fastq ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            [ meta_clone, fastq ]
    }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }


*/
