import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';

// General imports
import { HomeComponent } from './layout/home/home.component';
import { NotFoundComponent } from './meta/not-found/not-found.component';

// GOF app imports
import { GofHolderComponent } from './gof-holder/gof-holder.component';

//PhylogeneticsTrees Imports
import { PhylogeneticTreesComponent } from './phylogenetic-trees/phylogenetic-trees.component';
import {MsaViewerComponent} from './phylogenetic-trees/msa-viewer/msa-viewer.component';

const routes: Routes = [
  { path: 'home', component: HomeComponent },
  { path: 'loadGof', component: GofHolderComponent },
  { path: 'Trees', component: PhylogeneticTreesComponent },
  { path: 'MsaViewer', component: MsaViewerComponent },
  { path: '**', component: NotFoundComponent },
  { path: '', redirectTo: 'home', pathMatch: 'full' },

];

@NgModule({
  imports: [RouterModule.forRoot(routes)],
  exports: [RouterModule]
})
export class AppRoutingModule { }
