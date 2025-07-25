Publishing the data from your phantom calculations
==================================================================
Recommended best practice for open science is that parameter files, initial conditions
and snapshots from calculations with phantom should be made publicly available on publication.

FAIR Principles
----------------
According to the `FAIR principles for scientific data management <https://ardc.edu.au/resource/fair-data/>`__, your data should be:

- Findable, e.g. with links to and from the paper publishing the simulations
- Accessible, available for free in a publicly accessible repository
- Interoperable, data is labelled and able to be reused or converted
- Reusable, include enough information to be able to reproduce your simulations

Data curation
-------------
For calculations with phantom that have been published in a paper,
ideal practice is to upload the **entire calculation including .in and
.setup files, .ev files and all dump files in a public repository**.

See for example a dataset from Mentiplay et al. (2020) using figshare: `<https://doi.org/10.6084/m9.figshare.11595369.v1>`_

Or this example from Wurster, Bate & Price (2018) in the University of Exeter repository: `<https://doi.org/10.24378/exe.607>`_

However, size limitations may restrict preservation of all data, in which case we recommend saving:

- .in files
- .setup files
- .ev files
- dump files used to create figures in your paper, with a link to splash or sarracen in the metadata for how to read/convert these files
- dump files containing initial conditions, if these are non-trivial
- metadata including link to your publication or arXiv preprint, link to the phantom code, code version information and labelling of data corresponding to simulations listed in your paper

Zenodo community
----------------
To facilitate better data sharing between phantom users, we have set up a Zenodo community:

   https://zenodo.org/communities/phantom

Please join this community and let's learn from each other to create best-practice data curation. 
Zenodo currently has a 50Gb limit on data size, which is sufficient for the recommended list of files to save above.

Archiving your data to Google Drive using rclone
------------------------------------------------
You can use rclone to copy data from a remote cluster or supercomputing facility to Google Drive. This is not recommended as a long term storage solution but can facilitate short-term data sharing between users.

Set this up by logging into your supercomputer and typing::

   $ rclone config
   No remotes found - make a new one
   n) New remote
   s) Set configuration password
   q) Quit config

   name> dan-google-drive

   Storage> drive

   Google Application Client Id
   See https://rclone.org/drive/#making-your-own-client-id for how to create your own.

   client_id>   (enter client id you got from the instructions)
   client_secret>  (enter client secret you got from the instructions)

   scope> drive
   root_folder_id> (leave this blank)

   Edit advanced config? (y/n) n

   Use auto config?
   y/n> n

   Please go to the following link: https://accounts.google.com/o/oauth2/auth?access_type=offline&client_id=...
   (click on the link to approve)

   Configure this as a team drive?
   y/n> n

Check the above was successful by listing files on your remote drive using::

    $ rclone ls dan-google-drive:

To COPY files to your google drive, LEAVING a copy on the local machine, you can then use::

    $ rclone copy local_path remote_path

For example::

    $ rclone copy $HOME/runs/phantom/disc-test1 dan-google-drive:phantom/disc-test1
    
To MOVE files to your google drive and DELETE them from the cluster (e.g. to clear disc space)::

    $ rclone move $HOME/runs/phantom/disc-test1 dan-google-drive:phantom/disc-test1

To SYNC an entire directory tree onto your google drive, DELETING files ALSO ON THE DRIVE you can then use::

    $ rclone sync -i $HOME/runs dan-google-drive:runs

.. _authorship:

Authorship guidelines
======================
We can demonstrate that women are under-represented in co-authorships within the Phantom community, consistent with and inspired by `Joyce et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022PASP..134h4503J/abstract>`__. As Phantom is not a formal collaboration but rather a community of users, this disparity is likely due to unconscious bias. Because equality in our community is one of our core values, we seek to draw attention to this and change it.

We have drawn up the following guidelines for 'What constitutes a co-author' by surveying the participants at the 2025 European workshop. The responses have been grouped into 'do', 'might' and 'do not' constitute co-authorship from what the majority of the community voted. We recommend that lead authors use these guidelines to gauge whether to make someone a co-author:

**Does make co-authorship:**
- Analysed simulations or reduced observations
- Ran simulations or took observations
- Wrote some text for the paper
- Made a figure
- Provided proprietary observational data
- Had the original idea
- Wrote some code (that can't be cited/is not public)
- Supervised the lead author
- Worked hard on something that didn't end up working so wasn't included in the final paper

**Might make co-authorship:**
- Provided advice on the paper
- Provided advice on using code/analysis/observations
- Met regularly with you but didn't contribute new ideas
- Proof-read the paper and gave science suggestions

**Does not make co-authorship:**
- Provided funding
- Reciprocating co-authorship
- Proof-read the paper and gave grammar suggestions
- Part of your collaboration (but didn't do anything)
- Supervised a co-author (e.g. the student is a co-author and their supervisor wants to be a co-author too)
- When the paper was 95% done they gave you a teeny bit of advice
- They spoke to you briefly in a social setting about your idea
- They're your student but they didn't do anything constructive
- They are a senior person who lends credibility
- They're an expert on the topic but didn't do anything for the paper
- They are your friend
- Asked you a question after your talk at a conference

Other helpful information
--------------------------
- :doc:`General instructions for running on a remote cluster </getting-started/running-clusters>`
- `rclone userguide <https://rclone.org>`_
