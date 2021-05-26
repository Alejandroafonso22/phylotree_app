import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';

// Layout imports
import { HomeComponent } from './layout/home/home.component';
import { NotFoundComponent } from './meta/not-found/not-found.component';
import { FastaUploadComponent } from './layout/fasta-upload/fasta-upload.component';

// GOG app imports
import { GogHolderComponent } from './gog-holder/gog-holder.component';

const routes: Routes = [
  {path: 'home', component: HomeComponent},
  {path: 'loadGog', component: GogHolderComponent},
  {path: 'fasta', component: FastaUploadComponent},
  {path: '**', component: NotFoundComponent},
  {path: '', redirectTo: 'home', pathMatch: 'full'}
];

@NgModule({
  imports: [RouterModule.forRoot(routes)],
  exports: [RouterModule]
})
export class AppRoutingModule { }
