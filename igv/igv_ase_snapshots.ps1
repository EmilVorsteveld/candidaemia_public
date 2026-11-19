<#
IGV screenshots for the ASE analysis (GitHub issues #6 and #7): both Kinnex samples of a patient (D0, D5-7) and
the patient's lrWGS together in one view, as a gene overview and as a zoom on a variant.

Input: a targets CSV written by 08_allele_specific.Rmd, e.g. ase_example_genes_targets.csv, with columns
  name, view (gene/zoom), gene, individual, dna_id, kinnex_samples (";"-separated), locus (chr:start-end), snv (chr:pos)

The script finds each BAM on the T: drive (<root>\<id>\<date>\<id>.bam; with several dates it takes the newest),
writes one IGV batch file and, if -Igv is given, runs IGV on it. Without -Igv, open IGV and use
Tools > Run Batch Script on the batch file it prints.

Usage (PowerShell, from this folder):
  .\igv_ase_snapshots.ps1 -Targets .\ase_example_genes_targets.csv
  .\igv_ase_snapshots.ps1 -Targets .\ase_example_genes_targets.csv -Igv "C:\Program Files\IGV_2.19.5\igv.bat"
If PowerShell refuses to run scripts:  powershell -ExecutionPolicy Bypass -File .\igv_ase_snapshots.ps1 -Targets ...
#>
param(
  [Parameter(Mandatory = $true)][string]$Targets,
  [string]$KinnexRoot = "T:\PIgroup Alexander Hoischen\Emil\candidemia\kinnex\raw",
  [string]$WgsRoot = "T:\PIgroup Alexander Hoischen\Emil\candidemia\wgs\raw",
  [string]$OutDir = (Join-Path $PSScriptRoot "screenshots"),
  [string]$Genome = "hg38",
  [int]$PanelHeight = 350,   # max height (px) of each track panel in a screenshot; lower it if the view gets too tall
  [string]$Igv = ""          # path to igv.bat; empty = only write the batch file
)
$ErrorActionPreference = "Stop"

function Find-Bam([string]$Root, [string]$Id) {
  $dir = Join-Path $Root $Id
  $hits = @(Get-ChildItem -LiteralPath $dir -Filter "$Id.bam" -Recurse -File -ErrorAction SilentlyContinue | Sort-Object FullName -Descending)
  if ($hits.Count -eq 0) { Write-Warning "No BAM for $Id under $dir"; return $null }
  $bam = $hits[0].FullName
  $indexed = (Test-Path -LiteralPath "$bam.bai") -or (Test-Path -LiteralPath ($bam -replace '\.bam$', '.bai')) -or (Test-Path -LiteralPath "$bam.csi")
  if (-not $indexed) { Write-Warning "No index (.bai/.csi) for $bam; IGV needs an aligned, indexed BAM" }
  return $bam
}

$rows = Import-Csv -LiteralPath $Targets
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$batch = Join-Path $OutDir ("igv_batch_" + [IO.Path]::GetFileNameWithoutExtension($Targets) + ".txt")
$IgvPath = { param($p) '"' + ($p -replace '\\', '/') + '"' }   # IGV batch: forward slashes, quoted (paths contain spaces)

$lines = New-Object System.Collections.Generic.List[string]
$lines.Add("snapshotDirectory " + (& $IgvPath $OutDir))
$lines.Add("maxPanelHeight $PanelHeight")
$bamCache = @{}

# Load each patient's BAMs once, then take all of that patient's screenshots
foreach ($patient in ($rows | Group-Object individual)) {
  $first = $patient.Group[0]
  $ids = @($first.kinnex_samples -split ';' | Where-Object { $_ })
  $bams = @()
  foreach ($id in $ids) {
    if (-not $bamCache.ContainsKey($id)) { $bamCache[$id] = Find-Bam $KinnexRoot $id }
    if ($bamCache[$id]) { $bams += $bamCache[$id] }
  }
  $wgs = Find-Bam $WgsRoot $first.dna_id
  if ($wgs) { $bams += $wgs }
  if ($bams.Count -eq 0) { Write-Warning "Skipping $($patient.Name): no BAMs found"; continue }

  $lines.Add("new")
  $lines.Add("genome $Genome")
  foreach ($b in $bams) { $lines.Add("load " + (& $IgvPath $b)) }
  # Reads are grouped by haplotype tag (HP) where the BAM has one (haplotagged BAMs); otherwise they form one group
  $lines.Add("group TAG HP")

  foreach ($t in $patient.Group) {
    $lines.Add("goto $($t.locus)")
    if ($t.view -eq "zoom") {
      $chr, $pos = $t.snv -split ':'
      $lines.Add("region $chr $pos $pos")
      $lines.Add("sort base $($t.snv)")
      $lines.Add("squish")
    } else {
      $lines.Add("collapse")
    }
    $lines.Add("snapshot $($t.name).png")
  }
}
if ($Igv) { $lines.Add("exit") }

Set-Content -LiteralPath $batch -Value $lines -Encoding ASCII
Write-Host "Batch file: $batch"
Write-Host "Screenshots go to: $OutDir"
if ($Igv) {
  & $Igv -b $batch
} else {
  Write-Host "Open IGV and run it with Tools > Run Batch Script."
}
