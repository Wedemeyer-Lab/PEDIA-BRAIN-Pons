import numpy as np
import pandas as pd
import scipy.io
import scipy.stats
import scrublet as scr
import matplotlib.pyplot as plt
import os
import random
import sys

  
def scrub(sample_dir, input_dir, output_dir, gene_dir):
  if not os.path.isdir(output_dir): 
    os.mkdir(output_dir)
  # Set params for figures
  plt.rc('font', size=14)
  plt.rcParams['pdf.fonttype'] = 42
  
  #import sample names
  samples = pd.read_csv(sample_dir)
  samples = list(samples["sample_names"])
  

  # load counts matrix and genes
  for sample in samples:
    print("Running " + sample + "...")
    
    indir = input_dir + "/" + sample +"/outs"+ "/filtered_feature_bc_matrix/"
    counts_matrix = scipy.io.mmread(indir + "matrix.mtx.gz").T.tocsc()
    genes = np.array(scr.load_genes(gene_dir+sample+"_features.tsv", delimiter="\t", column=1))

    print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    print('Number of genes in gene list: {}'.format(len(genes)))

    # initialie Scrublet object
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

    # doublet_scores and predicted_doublets
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                              min_cells=3, 
                                                              min_gene_variability_pctl=85, 
                                                              n_prin_comps=30)
    scrub.call_doublets(threshold=0.25)    

    # plot doublet score histograms
    scrub.plot_histogram()
    plt.savefig(os.path.join(output_dir,sample+ '_doublet_score_histogram.png'))

    # get 2D embedding to visualize results
    print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))

    # # Uncomment to run tSNE - slow
    # print('Running tSNE...')
    # scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=0.9))

    # # Uncomment to run force layout - slow
    # print('Running ForceAtlas2...')
    # scrub.set_embedding('FA', scr.get_force_layout(scrub.manifold_obs_, n_neighbors=5. n_iter=1000))
          
    print('Done.')
      
    # plot doublet predictions on 2D embedding
    scrub.plot_embedding('UMAP', order_points=True);
    plt.savefig(os.path.join(output_dir, sample+'_doublet_score_umap.png'))

    # scrub.plot_embedding('tSNE', order_points=True);
    # scrub.plot_embedding('FA', order_points=True);

    # export doublet_scores, predicted_doublets as lists

    np.savetxt(os.path.join(output_dir, sample + '_doublet_scores.txt'), doublet_scores)

    np.savetxt(os.path.join(output_dir, sample + '_predicted_doublets.txt'), predicted_doublets)

         

if __name__ == "__main__":
  sample_dir = sys.argv[1]
  input_dir = sys.argv[2]
  output_dir = sys.argv[3]
  gene_dir = sys.argv[4]
  
  scrub(sample_dir=sample_dir, 
        input_dir=input_dir, 
        output_dir=output_dir, 
        gene_dir=gene_dir
        )

    
