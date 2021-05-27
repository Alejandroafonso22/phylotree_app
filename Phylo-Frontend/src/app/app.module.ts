import { NgModule } from '@angular/core';
import { BrowserModule } from '@angular/platform-browser';
import { NgxPaginationModule } from 'ngx-pagination'; // <-- import the module
import { AppRoutingModule } from './app-routing.module';
import { AppComponent } from './app.component';

// General imports
import { HomeComponent } from './layout/home/home.component';
import { FooterComponent } from './layout/footer/footer.component';
import { NavbarComponent } from './layout/navbar/navbar.component';
import { HeaderComponent } from './layout/header/header.component';
import { NotFoundComponent } from './meta/not-found/not-found.component';

// GOF app imports
import { GofHolderComponent } from './gof-holder/gof-holder.component';
import { MapViewComponent } from './gof-holder/map-view/map-view.component';

//PhyloGenetic-Trees imports 
import { PhylogeneticTreesComponent } from './phylogenetic-trees/phylogenetic-trees.component';
import { BackendConnectService } from './services/PhyloTrees/backend-connect.service';
import { HttpClientModule } from '@angular/common/http';
import { FormsModule, ReactiveFormsModule } from '@angular/forms';
import { AngularFileUploaderModule } from "angular-file-uploader";
import { TaxonFilterPipe } from './Pipes/taxon-filter.pipe';
import { MsaViewerComponent } from './phylogenetic-trees/msa-viewer/msa-viewer.component';


@NgModule({
  declarations: [
    AppComponent,
    GofHolderComponent,
    HomeComponent,
    FooterComponent,
    NavbarComponent,
    HeaderComponent,
    NotFoundComponent,
    PhylogeneticTreesComponent,
    TaxonFilterPipe,
    MsaViewerComponent,
  ],
  imports: [
    BrowserModule,
    AppRoutingModule,
    HttpClientModule,
    FormsModule,
    ReactiveFormsModule,
    NgxPaginationModule,
    AngularFileUploaderModule
  ],
  providers: [BackendConnectService],
  bootstrap: [AppComponent]
})
export class AppModule { }
