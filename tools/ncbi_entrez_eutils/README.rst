Galaxy NCBI Entrez Tools
========================

This repo requires a readme as administrators should very aware of some
restrictions NCBI places on the use of the Entrez service.

NCBI requests that you please limit large jobs to either weekends or
between 9:00 PM and 5:00 AM Eastern time during weekdays. This is not a
request that the Galaxy tool can easily service, so we've included it in
the disclaimer on every tool quite prominently.

Failure to comply with NCBI's policies may result in an block until
you/the user contacts NCBI and registers the tool ID and their email.

Note that these are *IP* level blocks so the Galaxy tools uses a
concatenation of the administrator's emails, and the user email, in
hopes that NCBI will contact all relevant parties should their system be
abused.

Additionally, since these are IP level blocks, the Galaxy tool author
(@erasche) recommends using the following ``jobs_conf.xml`` snippet in
order to place a system-wide restriction of 1 concurrent Entrez job
amongst all users.

.. code:: xml

    <destination id="entrez" runner="local">
    </destination>
    <limit type="concurrent_jobs" id="entrez">1</limit>
    <tools>
      <tool id="ncbi.eutils.efetch" destination="entrez" />
      <tool id="ncbi.eutils.esearch" destination="entrez" />
      <tool id="ncbi.eutils.epost" destination="entrez" />
      <tool id="ncbi.eutils.elink" destination="entrez" />
      <tool id="ncbi.eutils.einfo" destination="entrez" />
      <tool id="ncbi.eutils.esummary" destination="entrez" />
    </tools>

