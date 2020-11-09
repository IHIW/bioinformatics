We need to provide some guidelines for how vendors should provide reference sequences and coordinates, that we can enforce in submission of HML. This standard should allow clear interpretation but also should be clinically relevant if possible. This was the consensus from the discussion, a proposed set of rules for submission of HML for workshop projects. 

### Standard Validation Requirements

* Valid according to NMDP Gateway. 
* Valid to MIRING Standards
* Any GLStrings should be valid as well.

### Additional Requirements for Reference Sequences
          
* Known alleles can be reported against the known allele sequence in IMGT/HLA
* Allele Variants should be reported with a IMGT/HLA sequence as a reference
     * Reference will be against a full-length HLA allele sequence, including Intron Seuqences
     * Vendors are provided a short list of full-length alleles to use as reference sequences.
          * Selected to cover as wide of allele groups,
          * But should also represent structural variants (Insertions and Deletions)
          * Start with the list used in analysis from 17th IHIW (Kazu and Steve Mack list)
               * Alleles will be added to the list based on new sequences in IMGT/HLA
               * Based on sequence indels, or other needs
               * Or based on the availablity of full-length sequences.
     * If sequence is reported outside the reference, such as extended UTRs,
          * It's reported as Insertions at the beginning or the end of the reference.
* IMGT/HLA Database release version must be reported
     * should only on on a list of releases (3.40.0, 3.39.0, 3.38.0, etc)
     * At least for the last 5 (?) years.

### Potential Questions / Issues:
   
* Is this clinically Relevant for reporting novel allele? 
     * Is it more relevant to use a full-length common allele as a reference
          * Perhaps easier for reporting
     * Or the “closest allele" as a reference
          * Perhaps this is more relevant in the clinic
          * Ideally,  "variants" are reported against confirmed genotypes in glstring. This is what clinical labs are looking for. We would like to know where the variants are in the closest alleles. If the clinical techs find different variants in HML report, this will create a problem.

* How easy is it to translate back to IMGT/HLA for vendors?
     * Vendors often use internal databases with their own references.
     * It is work to re-align data to a required reference and report different coordinates

* Vendors may be hesitant to change sequence reporting
     * Especially if it messes with customer products
     * We may need to encourage workshop-specific, “research-only” tools for reporting

* How will vendors choose the "correct" reference sequence?
     * Can we provide simple tool for what allele to use?

* How will we make this easy to adapt?
     * We can make this proposal as a brief list of ideas, or
     * a complete set of example HML Files
          * Ranging from Simple to Complex
               * Examples should be based on the worst problems from the last workshop

* We need a system for vendors to check that the files are good.
     * Check sample files from vendors, both using validator software and by eye, so that we can understand the files.




