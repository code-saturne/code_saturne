<html>
<head>
  <meta http-equiv="content-type" content="text/html; charset=utf-8"/>
  <title></title>
</head>
<body

<big><b>Batch options for SLURM (Simple Linux Utility for Resource Management)</b></big>

<br><br>

For a more complete documentation, see SLURM documentation
 <a href="https://slurm.schedmd.com/documentation.html">(https://slurm.schedmd.com/documentation.html)</a>.<br>
The local cluster documentation may also include site-specific examples
and recommendations.<br>

<table>
  <thead>
    <tr>
      <th style="text-align:left">Option</th>
      <th style="text-align:left">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><samp>-t&lt;time&gt;, --time=&lt;time&gt;</samp></td>
      <td>walltime in minutes, minutes:seconds, hours:minutes:seconds, days-hours, days-hours:minutes, or days-hours:minutes:seconds</td>
    </tr>
    <tr>
      <td><samp>-N, --nodes=&lt;minnodes[-maxnodes]&gt;</samp></td>
      <td>number of allocated nodes</td>
    </tr>
    <tr>
      <td><samp>-n, --ntasks=&lt;ntasks&gt;</samp></td>
      <td>number of total tasks</td>
    </tr>
    <tr>
      <td><samp>--ntasks-per-node=&lt;ntasks&gt;</samp></td>
      <td>number of tasks per node</td>
    </tr>
    <tr>
      <td><samp>--ntasks-per-socket=&lt;ntasks&gt;</samp></td>
      <td>number of tasks per socket</td>
    </tr>
    <tr>
      <td><samp>--ntasks-per-core=&lt;ntasks&gt;</samp></td>
      <td>number of tasks per core</td>
    </tr>
    <tr>
      <td><samp>--cpus-per-task=&lt;nthreads&gt;</samp></td>
      <td>number of threads per task</td>
    </tr>
    <tr>
      <td><samp>--cpu-bind=cores, sockets</samp></td>
      <td>bind CPUs to cores or sockets</td>
    </tr>
    <tr>
      <td><samp>--mem-bind=local</samp></td>
      <td></td>
    </tr>
    <tr>
      <td><samp>--qos=&lt;qos&gt;</samp></td>
      <td>request given quality of service (for example "debug")</td>
    </tr>
    <tr>
      <td><samp>--contiguous</samp></th>
      <td>use contiguous nodes for minimal latency</td>
    </tr>
    <tr>
      <td><samp>--exclusive</samp></td>
      <td>do not share nodes with other jobs</td>
    </tr>
    <tr>
      <td><samp>--switches=&lt;count&gt;[@max-delay]</samp></td>
      <td>try to run on no more than count switches (for better performance)</td>
    </tr>
    <tr>
      <td><samp>--partition=&lt;name&gt;, -p&lt;name&gt;</samp></td>
      <td>partition (queue) (run "sinfo -s" for list of partitions)</td>
    </tr>
    <tr>
      <td><samp>--reservation=&lt;name&gt;</samp></td>
      <td>allocate resources from reservation</td>
    </tr>
    <tr>
      <td><samp>-o&lt;pattern&gt;, --output=&lt;pattern&gt;</samp></td>
      <td>output file name pattern</td>
    </tr>
    <tr>
      <td><samp>-e&lt;pattern&gt;, --error=&lt;pattern&gt;</samp></td>
      <td>error file name pattern</td>
    </tr>
    <tr>
      <td><samp>-J&lt;jobname&gt;, --job-name=&lt;jobname&gt;</samp></td>
      <td>job name</td>
    </tr>
  <tbody>
</table>

<br><br>

Lines containing SLURM options in a submitted script header should be prefixed by <samp>#SBATCH</samp>.
Multiple options may be used on a given line.<br>

<br>
<b>Example:</b>

<br><br>

<samp>
#SBATCH --nodes=2<br>
#SBATCH --ntasks-per-node=12<br>
#SBATCH --cpus-per-task=1<br>
#SBATCH --time=0:30:00<br>
#SBATCH --partition=par<br>
#SBATCH --output=job_%j.out.log --error=job_%j.err.log<br>
#SBATCH --job-name=nameandcase<br>
</samp>

</body>
</html>
