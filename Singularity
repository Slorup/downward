Bootstrap: docker
From: ubuntu:xenial

%post
    ## Install all necessary dependencies.
    apt-get update
    apt-get -y install cmake g++ make python autotools-dev automake gcc g++-multilib ca-certificates git
    rm -rf /var/lib/apt/lists/*

    git clone -b ipc-2018-seq-opt https://bitbucket.org/ipc2018-classical/team32.git /planner

    ## Build your planner
    cd /planner
    ./build.py release64 -j6

%runscript
    ## The runscript is called whenever the container is used to solve
    ## an instance.

    DOMAINFILE=$1
    PROBLEMFILE=$2
    PLANFILE=$3

    ## Call your planner.
    /planner/fast-downward.py \
        --build=release64 \
        --plan-file $PLANFILE \
        $DOMAINFILE \
        $PROBLEMFILE \
        --search "astar(cpdbs_symbolic(genetic_ss(use_ucb=true,num_episodes=10000000,num_collections=1,pdb_factory=symbolic,genetic_time_limit=900,time_limit=1.0,create_perimeter=true,use_first_goal_vars=true,use_norm_dist=true)))"

## Update the following fields with meta data about your submission.
## Please use the same field names and use only one line for each value.
%labels
Name        Complementary2
Description PDB heuristic which learns in situ how to adapt its pattern selection parameteres to maximize chance of finding "good" patterns as seach progresses and also complements any input heuristic.  This is an implementation of the complementary heuristic with an input perimeter PDB and symbolic PDBs, as described in the IJCAI paper: https://www.ijcai.org/proceedings/2017/601
Authors Santiago Franco<santiago.franco@gmail.com>, Levi H S Lelis<levilelis@gmail.com> and Mike W Barley<m.barley@auckland.ac.nz>
SupportsDerivedPredicates no
SupportsQuantifiedPreconditions no
SupportsQuantifiedEffects yes
