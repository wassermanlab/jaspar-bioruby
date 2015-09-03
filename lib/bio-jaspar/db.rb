#--
# = bio-jaspar/db.rb
# 
# Copyright:: (C) 2015-2015 Jessica Lee
# License:: Ruby License
# 
# JASPAR 2014 module
# 
# A direct import of Bio.motifs.jaspar.db module in Biopython
#++

require 'mysql2'

module Bio # :nodoc:
  module Jaspar
    # Class representing a JASPAR5 DB. The methods within are loosely based
    # on the perl TFBS::DB::JASPAR5 module.
    # 
    # This modules requires MySQLdb to be installed.
    #
    # *Note:* We will only implement reading of JASPAR motifs from the DB.
    # Unlike the perl module, we will not attempt to implement any methods to
    # store JASPAR motifs or create a new DB at this time.
    # 
    # <i>A direct import of Bio.motifs.jaspar module in Biopython</i>
    class JASPAR5

      # Default JASPAR collection to access
      JASPAR_DFLT_COLLECTION = 'CORE'

      # Unambiguous DNA bases
      DNA = Bio::Motifs::Alphabet.new.IUPAC_unambiguous_dna

      # Construct a JASPAR5 instance and connect to specified DB
      # 
      # === Arguments
      # +host+:: host name of the the JASPAR DB server
      # +name+:: name of the JASPAR database
      # +user+:: user name to connect to the JASPAR DB
      # +password+:: JASPAR DB password
      # 
      #--
      # TODO = Once BioRuby moves to support Ruby 2.x change this to 
      # keyword arguments
      #++
      def initialize(opts = {})
        opts = {
          :host => nil,
          :name => nil,
          :user => nil,
          :password => nil
        }.merge(opts)

        @name = opts[:name]
        @host = opts[:host]
        @user = opts[:user]
        @password = opts[:password]

        @dbh = Mysql2::Client.new(:host => opts[:host], 
                                  :username => opts[:user], 
                                  :password => opts[:password], 
                                  :database => opts[:name])
      end

      # Return a string represention of the JASPAR5 DB connection.
      # 
      # Format: <user>@<host>:<db_name>
      def to_s
        str = "#{@user}\@#{@host}:#{@name}"
        return str
      end

      # Fetch a single JASPAR motif from the DB by it's JASPAR matrix ID
      # (e.g. 'MA0001.1').
      # 
      # === Arguments
      # +id+:: JASPAR matrix ID. This may be a fully specified ID including the
      #        version number (e.g. MA0049.2) or just the base ID (e.g. MA0049).
      #        If only a base ID is provided, the latest version is returned.
      # 
      # === Returns
      # A Bio.motifs.Motif.jaspar object
      # 
      # *NOTE:* The perl TFBS module allows you to specify the type of matrix to
      # return (PFM, PWM, ICM) but matrices are always stored in JASAPR as
      # PFMs so this does not really belong here. Once a PFM is fetched the
      # pwm() and pssm() methods can be called to return the normalized and
      # log-odds matrices.
      def fetch_motif_by_id(id)
        # separate stable ID and version number
        base_id, version = Bio::Jaspar.split_jaspar_id(id)

        unless version
          # if ID contains no version portion, fetch the latest version
          version = _fetch_latest_version(base_id)
        end

        # fetch internal JASPAR matrix ID - also a check for validity
        int_id = nil
        if version
          int_id = _fetch_internal_id(base_id, version)
        end

        # fetch JASPAR motif using internal ID
        motif = nil
        if int_id
          motif = _fetch_motif_by_internal_id(int_id)
        end

        return motif
      end

      # Fetch a list of JASPAR motifs from a JASPAR DB by the given TF name(s).
      #
      # === Arguments
      # +name+:: a single name or list of names
      # 
      # === Returns
      # A list of Bio.motifs.Motif.japar objects
      # 
      # *NOTE:*
      # Names are not guaranteed to be unique. There may be more than one
      # motif with the same name. Therefore even if name specifies a single
      # name, a list of motifs is returned. This just calls
      # self.fetch_motifs(collection = None, tf_name = name).
      # 
      # This behaviour is different from the TFBS perl module's
      # get_matrix_by_name() method which always returns a single matrix,
      # issuing a warning message and returning the first matrix retrieved
      # in the case where multiple matrices have the same name.
      def fetch_motifs_by_name(name)
        return self.fetch_motifs(:collection => nil, :tf_name => name)
      end

      # Fetch a jaspar.Record (list) of motifs based on the provided selection
      # criteria.
      #
      # === Arguments
      # Except where obvious, all selection criteria arguments may be specified
      # as a single value or a list of values. Motifs must meet ALL the
      # specified selection criteria to be returned with the precedent
      # exceptions noted below.
      #
      # +all+::       Takes precedent of all other selection criteria.
      #               Every motif is returned. If 'all_versions' is also
      #               specified, all versions of every motif are returned,
      #               otherwise just the latest version of every motif is
      #               returned.
      # +matrix_id+:: Takes precedence over all other selection criteria except
      #               'all'.  Only motifs with the given JASPAR matrix ID(s)
      #               are returned. A matrix ID may be specified as just a base
      #               ID or full JASPAR IDs including version number. If only a
      #               base ID is provided for specific motif(s), then just the
      #               latest version of those motif(s) are returned unless
      #               'all_versions' is also specified.
      # +collection+:: Only motifs from the specified JASPAR collection(s)
      #                are returned. NOTE - if not specified, the collection
      #                defaults to CORE for all other selection criteria except
      #                'all' and 'matrix_id'. To apply the other selection
      #                criteria across all JASPAR collections, explicitly set
      #                collection=None.
      # +tf_name+::    Only motifs with the given name(s) are returned.
      # +tf_class+::   Only motifs of the given TF class(es) are returned.
      # +tf_family+::  Only motifs from the given TF families are returned.
      # +tax_group+:: Only motifs belonging to the given taxonomic supergroups
      #               are returned (e.g. 'vertebrates', 'insects', 'nematodes'
      #               etc.)
      # +species+::   Only motifs derived from the given species are returned.
      #               Species are specified as taxonomy IDs.
      # +data_type+:: Only motifs generated with the given data type (e.g.
      #               ('ChIP-seq', 'PBM', 'SELEX' etc.) are returned. NOTE -
      #               must match exactly as stored in the database.
      # +pazar_id+::  Only motifs with the given PAZAR TF ID are returned.
      # +medline+::   Only motifs with the given medline (PubmMed IDs) are 
      #               returned.
      # +min_ic+::    Only motifs whose profile matrices have at least this
      #               information content (specificty) are returned.
      # +min_length+:: Only motifs whose profiles are of at least this length
      #                are returned.
      # +min_sites+:: Only motifs compiled from at least these many binding
      #               sites are returned.
      # +all_versions+:: Unless specified, just the latest version of motifs
      #                  determined by the other selection criteria are returned
      #                  otherwise all versions of the selected motifs are
      #                  returned.
      #
      # === Returns
      # A Bio.motifs.Motif.jaspar.Record (list) of motifs.
      def fetch_motifs(opts = {})
        opts = {
          :collection => JASPAR_DFLT_COLLECTION, 
          :tf_name => nil, 
          :tf_class => nil, 
          :tf_family => nil, 
          :matrix_id => nil, 
          :tax_group => nil, 
          :species => nil, 
          :pazar_id => nil, 
          :data_type => nil, 
          :medline => nil, 
          :min_ic => 0, 
          :min_length => 0, 
          :min_sites => 0, 
          :all => false, 
          :all_versions => false
        }.merge(opts)

        # Fetch the internal IDs of the motifs using the criteria provided
        int_ids = _fetch_internal_id_list(
          :collection => opts[:collection],
          :tf_name => opts[:tf_name],
          :tf_class => opts[:tf_class],
          :tf_family => opts[:tf_family],
          :matrix_id => opts[:matrix_id],
          :tax_group => opts[:tax_group],
          :species => opts[:species],
          :pazar_id => opts[:pazar_id],
          :data_type => opts[:data_type],
          :medline => opts[:medline],
          :all => opts[:all],
          :all_versions => opts[:all_versions]
        )

        record = Bio::Jaspar::Record.new

        # Now further filter motifs returned above based on any specified
        # matrix specific criteria.
        int_ids.each do |int_id|
          motif = _fetch_motif_by_internal_id(int_id)

          # Filter motifs to those with matrix IC greater than min_ic
          if opts[:min_ic]
            if motif.pssm.mean < opts[:min_ic]
              next
            end
          end

          # Filter motifs to those with minimum length of min_length
          if opts[:min_length]
            if motif.length < opts[:min_length]
              next
            end
          end

          # Filter motifs to those composed of at least this many sites.
          # The perl TFBS module assumes column sums may be different but
          # this should be strictly enforced here we will ignore this and
          # just use the first column sum.
          if opts[:min_sites]
            num_sites = motif.alphabet.letters.map { |nt| motif.counts[nt][0] }.inject(:+)

            if num_sites < opts[:min_sites]
              next
            end
          end

          record << motif
        end

        return record
      end

      private

      # Get the latest version number for the given base_id
      def _fetch_latest_version(base_id)
        cur = @dbh.query("SELECT VERSION FROM MATRIX WHERE BASE_ID='#{base_id}'" \
                         " ORDER BY VERSION DESC LIMIT 1")
        row = cur.first

        latest = nil
        if row
          latest = row["VERSION"]
        else
          warn("Failed to fetch latest version number for JASPAR motif with" \
               " base ID #{base_id}. No JASPAR motif with this base ID" \
               " appears to exist in the databse.")
        end

        return latest
      end

      # Fetch the internal id for a base id + version. Also checks if this
      # combo exists or not
      def _fetch_internal_id(base_id, version)
        # fetch basic motif information
        cur = @dbh.query("SELECT ID FROM MATRIX WHERE BASE_ID='#{base_id}'" \
                         " and VERSION=#{version}")
        row = cur.first

        int_id = nil
        if row
          int_id = row["ID"]
        else
          warn("Failed to fetch internal database ID for JASPAR motif with" \
               " matrix ID '#{base_id}.#{version}'. No JASPAR motif with this" \
               " matrix ID appears to exist.")
        end

        return int_id
      end

      # Fetch basic motif information
      def _fetch_motif_by_internal_id(int_id)
        cur = @dbh.query("SELECT BASE_ID, VERSION, COLLECTION, NAME FROM" \
                         " MATRIX WHERE ID=#{int_id}")
        row = cur.first

        # This should never happen as it is an internal method. If it does
        # we should probably raise an exception
        unless row
          warn("Could not fetch JASPAR motif with internal ID = #{int_id}")
          return nil
        end

        base_id = row["BASE_ID"]
        version = row["VERSION"]
        collection = row["COLLECTION"]
        name = row["NAME"]

        matrix_id = [base_id, ".", version.to_s].join("")

        # fetch the counts matrix
        counts = _fetch_counts_matrix(int_id)

        # Create new JASPAR motif
        motif = Bio::Jaspar::Motif.new(matrix_id, 
                                       name,
                                       :counts => counts, 
                                       :collection => collection)

        # fetch species
        cur = @dbh.query("SELECT TAX_ID FROM MATRIX_SPECIES WHERE ID=#{int_id}")
        tax_ids = []
        cur.each { |r| tax_ids << r["TAX_ID"] }
        motif.species = tax_ids

        # fetch protein accession numbers
        cur = @dbh.query("SELECT ACC FROM MATRIX_PROTEIN WHERE ID=#{int_id}")
        accs = []
        cur.each { |r| accs << r["ACC"] }
        motif.acc = accs

        # fetch remaining annotation as tags from the ANNOTATION table
        cur = @dbh.query("SELECT TAG, VAL FROM MATRIX_ANNOTATION WHERE" \
                         " ID=#{int_id}")
        cur.each do |r|
          att = r["TAG"]
          val = r["VAL"]

          if att == "class"
            motif.tf_class = val
          elsif att == "family"
            motif.tf_family = val
          elsif att == "tax_group"
            motif.tax_group = val
          elsif att == "type"
            motif.data_type = val
          elsif att == "pazar_tf_id"
            motif.pazar_id = val
          elsif att == "medline"
            motif.medline = val
          elsif att == "comment"
            motif.comment = val
          else
            next
          end
        end

        return motif
      end

      # Fetch the counts matrix from the JASPAR DB by the internal ID
      #
      # Returns a Bio.motifs.matrix.GenericPositionMatrix
      def _fetch_counts_matrix(int_id)
        counts = {}

        DNA.letters.each do |base|
          base_counts = []
          cur = @dbh.query("SELECT VAL FROM MATRIX_DATA WHERE ID=#{int_id}" \
                           " and ROW='#{base}' ORDER BY COL")
          cur.each { |r| base_counts << r["VAL"] }
          counts[base] = base_counts.map(&:to_f)
        end

        return Bio::Motifs::GenericPositionMatrix.new(DNA, counts)
      end

      # Fetch a list of internal JASPAR motif IDs based on various passed
      # parameters which may then be used to fetch the rest of the motif data.
      #
      # Caller:
      #     fetch_motifs()
      #
      # Arguments:
      #     See arguments sections of fetch_motifs()
      #
      # Returns:
      #     A list of internal JASPAR motif IDs which match the given
      #     selection criteria arguments.
      #
      #
      # Build an SQL query based on the selection arguments provided.
      #
      # 1: First add table joins and sub-clauses for criteria corresponding to
      #    named fields from the MATRIX and MATRIX_SPECIES tables such as
      #    collection, matrix ID, name, species etc.
      #
      # 2: Then add joins/sub-clauses for tag/value parameters from the
      #    MATRIX_ANNOTATION table.
      #
      # For the surviving matrices, the responsibility to do matrix-based
      # feature filtering such as ic, number of sites etc, fall on the
      # calling fetch_motifs() method.
      def _fetch_internal_id_list(opts = {})
        opts = {
          :collection => JASPAR_DFLT_COLLECTION, 
          :tf_name => nil, 
          :tf_class => nil,
          :tf_family => nil, 
          :matrix_id => nil, 
          :tax_group => nil, 
          :species => nil,
          :pazar_id => nil, 
          :data_type => nil, 
          :medline => nil, 
          :all => false,
          :all_versions => false
        }.merge(opts)

        int_ids = []

        # Special case 1: fetch ALL motifs. Highest priority.
        # Ignore all other selection arguments.
        if opts[:all]
          cur = @dbh.query("SELECT ID FROM MATRIX")
          cur.each { |r| int_ids << r["ID"] }
          return int_ids
        end

        # Special case 2: fetch specific motifs by their JASPAR IDs. This
        # has higher priority than any other except the above 'all' case.
        # Ignore all other selection arguments.
        if opts[:matrix_id]
          # These might be either stable IDs or stable_ID.version.
          # If just stable ID and if all_versions == 1, return all versions,
          # otherwise just the latest
          if opts[:all_versions]
            opts[:matrix_id].each do |id|
              # ignore vesion here, this is a stupidity filter
              base_id, _ = Bio::Jaspar.split_jaspar_id(id)
              cur = @dbh.query("SELECT ID FROM MATRIX WHERE BASE_ID='#{base_id}'")
              cur.each{ |r| int_ids << r["ID"] }
            end
          else
            opts[:matrix_id].each do |id|
              # only the lastest version, or the requested version
              base_id, version = Bio::Jaspar.split_jaspar_id(id)

              unless version
                version = _fetch_latest_version(base_id)
              end

              int_id = nil
              if version
                int_id = _fetch_internal_id(base_id, version)
              end

              if int_id
                int_ids << int_id
              end
            end
          end

          return int_ids
        end

        tables = ["MATRIX AS m"]
        where_clauses = []

        # Select by MATRIX.COLLECTION
        if opts[:collection]
          if opts[:collection].is_a? Array
            # Multiple collections passed in as a list
            clause = "m.COLLECTION IN ('"
            clause += opts[:collection].join("','")
            clause += "')"
          else
            # A single collection - typical usage
            clause = "m.COLLECTION='#{opts[:collection]}'"
          end
          where_clauses << clause
        end

        # Select by MATRIX.NAME
        if opts[:tf_name]
          if opts[:tf_name].is_a? Array
            # Multiple names passed in as a list
            clause = "m.NAME IN ('"
            clause += opts[:tf_name].join("','")
            clause += "')"            
          else
            # A single name
            clause = "m.NAME='#{opts[:tf_name]}'"
          end
          where_clauses << clause
        end

        # Select by MATRIX_SPECIES.TAX_ID
        if opts[:species]
          tables << "MATRIX_SPECIES AS ms"
          where_clauses << "m.ID=ms.ID"

          # NOTE: species are numeric taxonomy IDs but stored as varchars
          # in the DB.
          if opts[:species].is_a? Array
            # Multiple tax IDs passed in as a list
            clause = "ms.TAX_ID IN ('"
            clause += opts[:species].join("','")
            clause += "')"            

          else
            # A single tax ID           
            clause = "m.TAX_ID='#{opts[:species]}'"
          end

          where_clauses << clause
        end

        # Tag based selection from MATRIX_ANNOTATION
        # Differs from perl TFBS module in that the matrix class explicitly
        # has a tag attribute corresponding to the tags in the database. This
        # provides tremendous flexibility in adding new tags to the DB and
        # being able to select based on those tags with out adding new code.
        # In the JASPAR Motif class we have elected to use specific attributes
        # for the most commonly used tags and here correspondingly only allow
        # selection on these attributes.

        # The attributes corresponding to the tags for which selection is
        # provided are:

        #    Attribute   Tag
        #    tf_class    class
        #    tf_family   family
        #    pazar_id    pazar_tf_id
        #    medline     medline
        #    data_type   type
        #    tax_group   tax_group
           
        # Select by TF class(es) (MATRIX_ANNOTATION.TAG="class")
        if opts[:tf_class]
          tables << "MATRIX_ANNOTATION AS ma1"
          where_clauses << "m.ID=ma1.ID"

          clause = "ma1.TAG='class'"
          if opts[:tf_class].is_a? Array
            # A list of TF classes
            clause += " AND ma1.VAL IN ('"
            clause += opts[:tf_class].join("','")
            clause += "')"
          else
            # A single TF class
            clause += " AND ma1.VAL='#{opts[:tf_class]}'"
          end

          where_clauses << clause
        end

        # Select by TF families (MATRIX_ANNOTATION.TAG="family")
        if opts[:tf_family]
          tables << "MATRIX_ANNOTATION AS ma2"
          where_clauses << "m.ID=ma2.ID"

          clause = "ma2.TAG='family'"
          if opts[:tf_family].is_a? Array
            # A list of TF families
            clause += " AND ma2.VAL IN ('"
            clause += opts[:tf_family].join("','")
            clause += "')"
          else
            # A single TF family
            clause += " AND ma2.VAL='#{opts[:tf_family]}'"
          end

          where_clauses << clause         
        end

        # Select by PAZAR TF ID(s) (MATRIX_ANNOTATION.TAG="pazar_tf_id")
        if opts[:pazar_id]
          tables << "MATRIX_ANNOTATION AS ma3"
          where_clauses << "m.ID=ma3.ID"

          clause = "ma3.TAG='pazar_tf_id'"
          if opts[:pazar_id].is_a? Array
            # A list of PAZAR IDs
            clause += " AND ma3.VAL IN ('"
            clause += opts[:pazar_id].join("','")
            clause += "')"
          else
            # A single PAZAR ID
            clause += " AND ma3.VAL='#{opts[:pazar_id]}'"
          end

          where_clauses << clause             
        end

        # Select by PubMed ID(s) (MATRIX_ANNOTATION.TAG="medline")
        if opts[:medline]
          tables << "MATRIX_ANNOTATION AS ma4"
          where_clauses << "m.ID=ma4.ID"

          clause = "ma4.TAG='medline'"
          if opts[:medline].is_a? Array
            # A list of PubMed IDs
            clause += " AND ma4.VAL IN ('"
            clause += opts[:medline].join("','")
            clause += "')"
          else
            # A single PubMed ID
            clause += " AND ma4.VAL='#{opts[:medline]}'"
          end

          where_clauses << clause                 
        end

        # Select by data type(s) used to compile the matrix
        # (MATRIX_ANNOTATION.TAG="type")
        if opts[:data_type]
          tables << "MATRIX_ANNOTATION AS ma5"
          where_clauses << "m.ID=ma5.ID"

          clause = "ma5.TAG='type'"
          if opts[:data_type].is_a? Array
            # A list of data types
            clause += " AND ma5.VAL IN ('"
            clause += opts[:data_type].join("','")
            clause += "')"
          else
            # A single data type
            clause += " AND ma5.VAL='#{opts[:data_type]}'"
          end

          where_clauses << clause           
        end

        # Select by taxonomic supergroup(s) (MATRIX_ANNOTATION.TAG="tax_group")
        if opts[:tax_group]
          tables << "MATRIX_ANNOTATION AS ma6"
          where_clauses << "m.ID=ma6.ID"

          clause = "ma6.TAG='tax_group'"
          if opts[:tax_group].is_a? Array
            # A list of tax IDs
            clause += " AND ma6.VAL IN ('"
            clause += opts[:tax_group].join("','")
            clause += "')"
          else
            # A single tax ID
            clause += " AND ma6.VAL='#{opts[:tax_group]}'"
          end

          where_clauses << clause
        end

        sql = "SELECT DISTINCT(m.ID) FROM " + tables.join(", ")

        if where_clauses
          sql += " WHERE " + where_clauses.join(" AND ")
        end

        cur = @dbh.query(sql)

        cur.each do |r|
          id = r["ID"]
          if opts[:all_versions]
            int_ids << id
          else
            # is the latest version?
            if _is_latest_version(id)
              int_ids << id
            end
          end
        end

        if int_ids.length < 1
          warn("Zero motifs returned with current select criteria")
        end

        return int_ids
      end

      # Does this internal ID represent the latest version of the JASPAR
      # matrix (collapse on base ids)     
      def _is_latest_version(int_id)
        cur = @dbh.query("SELECT COUNT(*) FROM MATRIX WHERE" \
                         " BASE_ID=(SELECT BASE_ID FROM MATRIX WHERE" \
                         " ID=#{int_id}) AND VERSION>(SELECT VERSION FROM" \
                         " MATRIX WHERE ID=#{int_id})")
        
        row = cur.first
        
        count = row["COUNT(*)"]
        
        if count == 0
          # no matrices with higher version ID and same base id
          return true
        end

        return false
      end

    end
  end
end